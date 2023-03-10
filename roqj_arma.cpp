#include "roqj.h"
// ------------------------- FUNCTIONS DEFINITIONS -------------------------
bool isNormalised (const cx_vec &psi) {return arma::norm(psi) == 1;}

cx_mat comm (const cx_mat &A, const cx_mat &B) {return A*B-B*A;}

cx_mat anticomm (const cx_mat &A, const cx_mat &B) {return A*B+B*A;}

cx_mat projector (const cx_vec &psi) {return psi*psi.t();}

cx_mat BlochToMatrix (double x, double y, double z) {
  double r = sqrt(x*x + y*y + z*z);
  if (r > 1.) {x /= r; y /= r; z /= r;}
  complex<double> I(0,1), one(1,0);
  cx_mat sigma_z = {{one,0},{0,-one}}, Id = {{one,0},{0,one}}, sigma_x = {{0,one},{one,0}}, sigma_y = {{0,-I},{I,0}};
  return .5*(Id + x*sigma_x + y*sigma_y + z*sigma_z);
}

// ------------------------- METHODS DEFINITIONS -------------------------
// ------------------------- ROQJ class -------------------------
// --- Constructors
void roqj::initialize (int N_ensemble, double t_i, double t_f, double dt, int N_copies, int dim_Hilbert_space, bool print_trajectory, int N_traj_print, bool verbose, double threshold) {
  set_N_ensemble(N_ensemble);
  set_time(t_i, t_f, dt);
  set_N_copies(N_copies);
  set_dim_Hilbert_space(dim_Hilbert_space);
  _print_trajectory = print_trajectory;
  set_N_traj_print(N_traj_print);
  set_initial_state();
  _verbose=verbose;
  set_threshold(threshold);
}

roqj::roqj (int N_ensemble, double t_i, double t_f, double dt, int N_copies, int dim_Hilbert_space, bool print_trajectory, int N_traj_print, bool verbose, double threshold) {
  initialize(N_ensemble, t_i, t_f, dt, N_copies, dim_Hilbert_space, print_trajectory, N_traj_print,verbose,threshold);
  srand(time(NULL));
}

// --- Setter
void roqj::set_N_ensemble (int N_ensemble) {
  if (N_ensemble <= 0) _N_ensemble = N_ensemble_default;
  else _N_ensemble = N_ensemble;
}

void roqj::set_N_copies (int N_copies) {
  if (N_copies <= 0) _N_copies = N_copies_default;
  else _N_copies = N_copies;
}

void roqj::set_t_i (double t_i) {
  if (t_i >= _t_f || t_i < 0.) _t_i = t_f_default;
  else _t_i = t_i;
}

void roqj::set_t_f (double t_f) {
  if (t_f <= _t_i) {
    if (t_f_default > _t_i) _t_f = t_f_default;
    else _t_f = 10.*_t_i;
  }
  else _t_f = t_f;
}

void roqj::set_dt (double dt) {
  if (dt == dt_default || dt <= 0. || dt >= _t_f - _t_i) _dt = _t_f/10000.;
  else _dt = dt;
  _num_timesteps = (int)(_t_f - _t_i)/_dt;
  _observable.resize(_num_timesteps);
  _sigma_observable.resize(_num_timesteps);
}

void roqj::set_time (double t_i, double t_f, double dt) {
  if (t_f > t_i) {
    _t_i = t_i;
    _t_f = t_f;
  }
  else if (t_f == t_i) {
    _t_i = t_i_default;
    _t_f = t_f_default;
  }
  else {
    _t_i = t_f;
    _t_f = t_i;
  }
  set_dt(dt);
}

void roqj::set_dim_Hilbert_space (int dim_Hilbert_space) {
  if (dim_Hilbert_space < 2) _dim_Hilbert_space = dim_Hilbert_space_default;
  else _dim_Hilbert_space = dim_Hilbert_space;
}

void roqj::set_N_traj_print (int N_traj_print) {
  if (_print_trajectory && (N_traj_print <= 0 || N_traj_print > _N_ensemble)) _N_traj_print = N_traj_print_default;
  else if (!_print_trajectory) _N_traj_print = 0;
  else _N_traj_print = N_traj_print;
}

void roqj::set_initial_state (const cx_vec &psi_i) {
  double n = arma::norm(psi_i);
  if (n == 0 || psi_i.size() != _dim_Hilbert_space) {
    set_initial_state();
    return;
  }
  _initial_state = normalise(psi_i);
}

void roqj::set_initial_state () {
  _initial_state = normalise(cx_vec(_dim_Hilbert_space, arma::fill::ones));
}

void roqj::set_print_traj (bool print) {_print_trajectory = print;}

void roqj::set_verbose (bool verbose) {_verbose = verbose;}

void roqj::set_threshold (double threshold) {
  if (threshold < 0.)
    threshold = -threshold;
  if (threshold > 0.1)
    _threshold = threshold_default;
  else _threshold = threshold;
}

// --- Getters
int roqj::get_N_ensemble () const {return _N_ensemble;}
int roqj::get_N_copies () const {return _N_copies;}
int roqj::get_dim_Hilbert_space () const {return _dim_Hilbert_space;}
int roqj::get_N_traj_print () const {return _N_traj_print;}
double roqj::get_t_i () const {return _t_i;}
double roqj::get_t_f () const {return _t_f;}
double roqj::get_dt () const {return _dt;}
double roqj::get_threshold () const {return _threshold;}
cx_vec roqj::get_initial_state () const {return _initial_state;}
vec roqj::get_observable () const {return _observable;}
vec roqj::get_error_observable () const {return _sigma_observable;}

vec roqj::get_observable (string file_out) const {
  ofstream out;
  out.open(file_out);
  for (int i = 0; i < _num_timesteps; ++i)
    out << _observable[i] << endl;
  return _observable;
}

vec roqj::get_error_observable (string file_out) const {
  ofstream out;
  out.open(file_out);
  for (int i = 0; i < _num_timesteps; ++i)
    out << _sigma_observable[i] << endl;
  return _sigma_observable;
}


// --- Run single iteration
vec roqj::run_single_iterations (bool verbose) const {
  vec observables(_num_timesteps);
  int n_observable = 0;

  // Allocating _N_ensemble copies of the initial state
  std::vector<cx_vec> psi(_N_ensemble);
  for (int i = 0; i <= _N_ensemble; ++i)
    psi[i] = _initial_state;

  // Exact solution
  cx_mat rho_ex(_dim_Hilbert_space, _dim_Hilbert_space);
  ofstream out_ex, traj;
  if (verbose) {
    rho_ex = projector(_initial_state);
    out_ex.open("analytic.txt");
    traj.open("trajectories.txt");
  }
  
  
  // Time evolution
  for (double t = _t_i; t <= _t_f; t += _dt) {
    // Prints and evolves the exact solution
    if (verbose) {
      out_ex << observable(rho_ex) << endl;
      rho_ex = rho_ex + (-complex<double>(0,1)*comm(H(t),rho_ex) + J(rho_ex,t) - 0.5*anticomm(Gamma(t),rho_ex))*_dt;
      rho_ex /= arma::trace(rho_ex);
    }

    // Average state
    cx_mat rho(_dim_Hilbert_space, _dim_Hilbert_space, arma::fill::zeros);

    // Cycle on the ensemble members
    for (int i = 0; i < _N_ensemble; ++i) {
      // Prints the trajectories
      if (verbose && i < _N_traj_print && _print_trajectory)
        traj << observable(projector(psi[i])) << " ";

      // Updates the average
      rho += projector(psi[i])/((double)_N_ensemble);

      cx_mat R = J(projector(psi[i]),t) + 0.5*(C(projector(psi[i]), t)*projector(psi[i]) + projector(psi[i])*C(projector(psi[i]), t).t());
      
      // Draws a random number and calculates whether the evolution is deterministic or via a jump
      double z = (double)rand()/((double)RAND_MAX);

      if (z < real(arma::trace(R))*_dt) // Jump
        psi[i] = this->jump(R,z);
      else {// Free evolution
        cx_mat K = H(t) + 0.5*(imag(C(projector(psi[i]), t)) - complex<double>(0.,1.)*(Gamma(t) + real(C(projector(psi[i]), t)) ) );
        psi[i] -= K*psi[i]*complex<double>(0.,1.)*_dt;
      }
      psi[i] = normalise(psi[i]);
    }
    // Storing the observable
    observables[n_observable] = observable(rho);
    n_observable++;

    if (verbose) traj << endl;
  }
  return observables;
}


// --- Running with all the copies
void roqj::run () {
  if (_verbose)
    print_info();

  ofstream params;
  params.open("params.txt");
  params << _N_copies << endl << _N_ensemble << endl << _t_i << endl << _t_f << endl << _dt << endl << _print_trajectory << endl << _N_traj_print << endl << _dim_Hilbert_space;
  params.close();

  arma::mat matrix_observables(_num_timesteps, _N_copies, arma::fill::zeros);
  for (int i = 0; i < _N_copies; ++i) {
    if (_verbose)
      cout << "Running copy " << i+1 << "/" << _N_copies << "...\n";
    vec this_obs = run_single_iterations(i==0);
    for (int j = 0; j < _num_timesteps; ++j)
      matrix_observables(j,i) = this_obs(j);
  }

  for (int i = 0; i < _num_timesteps; ++i) {
    _observable[i] = mean(matrix_observables.row(i));
    _sigma_observable[i] = stddev(matrix_observables.row(i));
  }
}

cx_vec roqj::jump (const cx_mat &R, double z) const {
  cx_vec eigval(_dim_Hilbert_space);
  cx_mat eigvec(_dim_Hilbert_space,_dim_Hilbert_space);
  eig_gen(eigval, eigvec, R);

  // Chose in which eigenvalue perform the jump
  double sum_previous_eigs = 0.;
  for (int j = 0; j < _dim_Hilbert_space; ++j) {
    if (real(eigval[j]) < -_threshold) {
      cerr << "Negative rate - reverse jump. NOT IMPLEMENTED\n";
      exit(EXIT_FAILURE);
    }
    // If z is in the j-th bin, it jumps to the j-th eigenstate
    if (z >= sum_previous_eigs*_dt && z < (sum_previous_eigs + real(eigval[j]))*_dt)
      return eigvec.col(j)*exp(-arg(eigvec.col(j)[1]));
    sum_previous_eigs += real(eigval[j]);
  }
  return cx_vec(_dim_Hilbert_space,arma::fill::ones);
}

void roqj::reset () {
  _observable.reset();
  _sigma_observable.reset();
}

void roqj::print_info () const {
  cout << "\nRate Operator Quantum Jumps - running " << _N_copies << " copies.\n";
  cout << "\tEnsemble size = " << _N_ensemble << ", " << _dim_Hilbert_space << "-dimensional Hilbert space,\n";
  cout << "\tt_i = " << _t_i << ", t_f = " << _t_f << ", dt = " << _dt << ",\n";
  if (_print_trajectory)
    cout << "\tPrinting " << _N_traj_print << " trajectories.\n\n";
  else cout << endl;
}











// ------------------------- ROQJ class -------------------------
// --- Constructors
qubit_roqj::qubit_roqj (int N_ensemble, double t_i, double t_f, double dt, int N_copies, bool print_trajectory, int N_traj_print, bool verbose) {
  srand(time(NULL));
  initialize(N_ensemble, t_i, t_f, dt, N_copies, 2, print_trajectory, N_traj_print, verbose);
}

// -- Set initial state vector
void qubit_roqj::set_initial_state (const cx_vec &psi) {
  if (arma::norm(psi) == 0 || psi.size() != 2) {
    set_initial_state();
    return;
  }
  _initial_state = normalise(psi);
}

// Default initial state
void qubit_roqj::set_initial_state () {_initial_state = {\./sqrt(2.), 1./sqrt(2.)};}

// --- Jump
cx_vec qubit_roqj::jump (const cx_mat &R, double z) const {
  cx_vec eigval(2);
  cx_mat eigvec(2,2);
  eig_gen(eigval, eigvec, R);

  double lambda1 = real(eigval[0]), lambda2 = real(eigval[1]), pjump1 = lambda1*_dt;
  if (lambda1 >= -_threshold && lambda2 >= -_threshold) {// Normal jump
    // With probability pjump1, it jumps to the first eigenstate of R
    if (z <= pjump1)
      return eigvec.col(0);
    else return eigvec.col(1);
  }
  else {// Reverse jump ----- Not implemented??
    cerr << "Negative rate - reverse jump. NOT IMPLEMENTED\n";
    exit(EXIT_FAILURE);
  }
  return cx_vec(2,arma::fill::ones);
}


// --- Run single iteration
vec qubit_roqj::run_single_iterations (bool verbose) const {
  vec observables(_num_timesteps);
  int n_observable = 0;

  // Allocating _N_ensemble copies of the initial state
  std::vector<cx_vec> psi(_N_ensemble);
  for (int i = 0; i <= _N_ensemble; ++i)
    psi[i] = _initial_state;

  // Exact solution
  cx_mat rho_ex(2, 2);
  ofstream out_ex, traj;
  if (verbose) {
    rho_ex = projector(_initial_state);
    out_ex.open("analytic.txt");
    traj.open("trajectories.txt");
  }
  
  
  // Time evolution
  for (double t = _t_i; t <= _t_f; t += _dt) {
    // Prints and evolves the exact solution
    if (verbose) {
      out_ex << observable(rho_ex) << endl;
      rho_ex = rho_ex + (-complex<double>(0,1)*comm(H(t),rho_ex) + J(rho_ex,t) - 0.5*anticomm(Gamma(t),rho_ex))*_dt;
      rho_ex /= arma::trace(rho_ex);
    }

    // Average state
    cx_mat rho(2, 2, arma::fill::zeros);

    // Cycle on the ensemble members
    for (int i = 0; i < _N_ensemble; ++i) {
      // Prints the trajectories
      if (verbose && i < _N_traj_print && _print_trajectory)
        traj << observable(projector(psi[i])) << " ";

      // Updates the average
      rho += projector(psi[i])/((double)_N_ensemble);

      cx_mat R = J(projector(psi[i]),t) + 0.5*(C(projector(psi[i]), t)*projector(psi[i]) + projector(psi[i])*C(projector(psi[i]), t).t());
      
      // Draws a random number and calculates whether the evolution is deterministic or via a jump
      double z = (double)rand()/((double)RAND_MAX);

      if (z < real(arma::trace(R))*_dt) // Jump
        psi[i] = this->jump(R,z);
      else {// Free evolution
        cx_mat K = H(t) + 0.5*(imag(C(projector(psi[i]), t)) - complex<double>(0.,1.)*(Gamma(t) + real(C(projector(psi[i]), t)) ) );
        psi[i] -= K*psi[i]*complex<double>(0.,1.)*_dt;
      }
      psi[i] = normalise(psi[i]);
    }
    // Storing the observable
    observables[n_observable] = observable(rho);
    n_observable++;

    if (verbose) traj << endl;
  }
  return observables;
}


// --- Running with all the copies
void qubit_roqj::run () {
  if (_verbose)
    print_info();

  ofstream params;
  params.open("params.txt");
  params << _N_copies << endl << _N_ensemble << endl << _t_i << endl << _t_f << endl << _dt << endl << _print_trajectory << endl << _N_traj_print << endl << 2;
  params.close();

  arma::mat matrix_observables(_num_timesteps, _N_copies, arma::fill::zeros);
  for (int i = 0; i < _N_copies; ++i) {
    if (_verbose)
      cout << "Running copy " << i+1 << "/" << _N_copies << "...\n";
    vec this_obs = run_single_iterations(i==0);
    for (int j = 0; j < _num_timesteps; ++j)
      matrix_observables(j,i) = this_obs(j);
  }

  for (int i = 0; i < _num_timesteps; ++i) {
    _observable[i] = mean(matrix_observables.row(i));
    _sigma_observable[i] = stddev(matrix_observables.row(i));
  }
}
