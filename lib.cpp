#include "lib.h"

// Free-evolution efective Hamiltonian
cx_mat K (const cx_mat &rho, double t) {
  return H(t) + 0.5*(imag(C(rho, t)) - complex<double>(0.,1.)*(Gamma(t) + real(C(rho, t)) ) );
}

cx_mat J_prime (const cx_mat &rho, double t) {
  return J(rho,t) + 0.5*(C(rho, t)*rho + rho*C(rho, t).t());
}

// ------------------------- FUNCTIONS DEFINITIONS -------------------------
bool isNormalised (const cx_vec &psi) {return arma::norm(psi) == 1;}

cx_mat comm (const cx_mat &A, const cx_mat &B) {return A*B-B*A;}

cx_mat anticomm (const cx_mat &A, const cx_mat &B) {return A*B+B*A;}

cx_mat projector (const cx_vec &psi) {return psi*psi.t();}

// ------------------------- METHODS DEFINITIONS -------------------------
// --- Constructors
void roqj::initialize (int N_ensemble, double t_i, double t_f, double dt, int N_copies, int dim_Hilbert_space, bool print_trajectory, int N_traj_print) {
  set_N_ensemble(N_ensemble);
  set_time(t_i, t_f, dt);
  set_N_copies(N_copies);
  set_dim_Hilbert_space(dim_Hilbert_space);
  _print_trajectory = print_trajectory;
  set_N_traj_print (N_traj_print);
  set_initial_state();
}

roqj::roqj (int N_ensemble, double t_i, double t_f, double dt, int N_copies, int dim_Hilbert_space, bool print_trajectory, int N_traj_print) {
  initialize(N_ensemble, t_i, t_f, dt, N_copies, dim_Hilbert_space, print_trajectory, N_traj_print);
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
  if (dt == dt_default || dt < 0. || dt >= _t_f - _t_i) _dt = dt = _t_f/10000.;
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

// --- Getters
int roqj::get_N_ensemble () const {return _N_ensemble;}
int roqj::get_N_copies () const {return _N_copies;}
int roqj::get_dim_Hilbert_space () const {return _dim_Hilbert_space;}
int roqj::get_N_traj_print () const {return _N_traj_print;}
double roqj::get_t_i () const {return _t_i;}
double roqj::get_t_f () const {return _t_f;}
double roqj::get_dt () const {return _dt;}
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
      rho_ex += (-complex<double>(0,1)*comm(H(t),rho_ex) + J(rho_ex,t) - 0.5*anticomm(Gamma(t),rho_ex))*_dt;
    }

    // Average state
    cx_mat rho(_dim_Hilbert_space, _dim_Hilbert_space, arma::fill::zeros);

    // Cycle on the ensemble members
    for (int i = 0; i < _N_ensemble; ++i) {
      // Prints the trajectories
      if (verbose && i < _N_traj_print)
        traj << observable(projector(psi[i])) << " ";

      // Updates the average
      rho += projector(psi[i])/((double)_N_ensemble);

      cx_mat R = J_prime(projector(psi[i]),t);
      
      // Draws a random number and calculates whether the evolution is deterministic or via a jump
      double z = (double)rand()/((double)RAND_MAX);

      if (z < real(arma::trace(R))*_dt) { // Jump
        cx_vec eigval(_dim_Hilbert_space);
        cx_mat eigvec(_dim_Hilbert_space,_dim_Hilbert_space);
        eig_gen(eigval, eigvec, R);

        // Chose in which eigenvalue perform the jump
        double sum_previous_eigs = 0.;
        bool already_jumped = false;
        for (int j = 0; j < _dim_Hilbert_space && !already_jumped; ++j) {
          // If z is in the j-th bin, it jumps to the j-th eigenstate
          if (z >= sum_previous_eigs*_dt && z < (sum_previous_eigs + real(eigval[j]))*_dt) {
            already_jumped = true;
            psi[i] = eigvec.col(j);
          }
          sum_previous_eigs += real(eigval[j]);
        }
      }
      else { // Free evolution
        psi[i] -= K(rho, t)*psi[i]*complex<double>(0.,1.)*_dt;
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
  cout << "\nRate Operator Quantum Jumps - running " << _N_copies << " copies.\n";
  cout << "\tEnsemble size = " << _N_ensemble << ", " << _dim_Hilbert_space << "-dimensional Hilbert space,\n";
  cout << "\tt_i = " << _t_i << ", t_f = " << _t_f << ", dt = " << _dt << ",\n";
  if (_print_trajectory)
    cout << "\tPrinting " << _N_traj_print << " trajectories.\n\n";
  else cout << endl;

  ofstream params;
  params.open("params.txt");
  params << _N_copies << endl << _N_ensemble << endl << _t_i << endl << _t_f << endl << _dt << endl << _print_trajectory << endl << _N_traj_print << endl << _dim_Hilbert_space;
  params.close();

  arma::mat matrix_observables(_num_timesteps, _N_copies, fill::zeros);
  for (int i = 0; i < _N_copies; ++i) {
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

void roqj::reset () {
  _observable.reset();
  _sigma_observable.reset();
}