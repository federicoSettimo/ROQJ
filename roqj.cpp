#include "roqj.h"
// ------------------------- FUNCTIONS DEFINITIONS -------------------------
bool isNormalized (const VectorXcd &psi) {return psi.norm() == 1;}

MatrixXcd comm (const MatrixXcd &A, const MatrixXcd &B) {return A*B-B*A;}

MatrixXcd anticomm (const MatrixXcd &A, const MatrixXcd &B) {return A*B+B*A;}

MatrixXcd projector (const VectorXcd &psi) {return psi*psi.transpose().conjugate();}

MatrixXcd BlochToMatrix (double x, double y, double z) {
  double r = sqrt(x*x + y*y + z*z);
  if (r > 1.) {x /= r; y /= r; z /= r;}
  return .5*(id + x*sigma_x + y*sigma_y + z*sigma_z);
}

// ------------------------- METHODS DEFINITIONS -------------------------
// ------------------------- ROQJ class -------------------------
// --- Constructors
void roqj::initialize (int N_states, double t_i, double t_f, double dt, int N_copies, int dim_Hilbert_space, bool print_trajectory, int N_traj_print, bool verbose, double threshold) {
  set_N_states(N_states);
  set_time(t_i, t_f, dt);
  set_N_copies(N_copies);
  set_dim_Hilbert_space(dim_Hilbert_space);
  _print_trajectory = print_trajectory;
  set_N_traj_print(N_traj_print);
  set_initial_state();
  _verbose=verbose;
  set_threshold(threshold);
}

roqj::roqj (int N_states, double t_i, double t_f, double dt, int N_copies, int dim_Hilbert_space, bool print_trajectory, int N_traj_print, bool verbose, double threshold) {
  initialize(N_states, t_i, t_f, dt, N_copies, dim_Hilbert_space, print_trajectory, N_traj_print,verbose,threshold);
  srand(0);
}

// --- Setter
void roqj::set_N_states (int N_states) {
  if (N_states <= 0) _N_states = N_states_default;
  else _N_states = N_states;
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
  _num_timesteps = (int)(_t_f - _t_i)/_dt + 1;
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
  if (_print_trajectory && (N_traj_print <= 0 || N_traj_print > _N_states)) _N_traj_print = N_traj_print_default;
  else if (!_print_trajectory) _N_traj_print = 0;
  else _N_traj_print = N_traj_print;
}

void roqj::set_initial_state (const VectorXcd &psi_i) {
  double n = psi_i.norm();
  if (n == 0 || psi_i.size() != _dim_Hilbert_space) {
    set_initial_state();
    return;
  }
  _initial_state = psi_i.normalized();
}

void roqj::set_initial_state () {
  _initial_state = VectorXcd::Ones(_dim_Hilbert_space).normalized();
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
int roqj::get_N_states () const {return _N_states;}
int roqj::get_N_copies () const {return _N_copies;}
int roqj::get_dim_Hilbert_space () const {return _dim_Hilbert_space;}
int roqj::get_N_traj_print () const {return _N_traj_print;}
double roqj::get_t_i () const {return _t_i;}
double roqj::get_t_f () const {return _t_f;}
double roqj::get_dt () const {return _dt;}
double roqj::get_threshold () const {return _threshold;}
VectorXcd roqj::get_initial_state () const {return _initial_state;}
VectorXd roqj::get_observable () const {return _observable;}
VectorXd roqj::get_error_observable () const {return _sigma_observable;}

VectorXd roqj::get_observable (string file_out) const {
  ofstream out;
  out.open(file_out);
  for (int i = 0; i < _num_timesteps; ++i)
    out << _observable[i] << endl;
  return _observable;
}

VectorXd roqj::get_error_observable (string file_out) const {
  ofstream out;
  out.open(file_out);
  for (int i = 0; i < _num_timesteps; ++i)
    out << _sigma_observable[i] << endl;
  return _sigma_observable;
}


// --- Run single iteration
VectorXd roqj::run_single_iterations (bool verbose) const {
  VectorXd observables(_num_timesteps);
  int n_observable = 0;

  // Allocating _N_states copies of the initial state
  std::vector<VectorXcd> psi;
  for (int i = 0; i <= _N_states; ++i)
    psi.push_back(_initial_state);

  // Exact solution
  MatrixXcd rho_ex(_dim_Hilbert_space, _dim_Hilbert_space);
  rho_ex = projector(_initial_state);
  ofstream out_ex, traj;
  if (verbose) {
    out_ex.open("analytic.txt");
    traj.open("trajectories.txt");
  }
  
  
  // Time evolution
  for (double t = _t_i; t <= _t_f; t += _dt) {
    // Prints and evolves the exact solution
    if (verbose) {
      out_ex << observable(rho_ex) << endl;
      rho_ex = rho_ex + (-complex<double>(0,1)*comm(H(t),rho_ex) + J(rho_ex,t) - 0.5*anticomm(Gamma(t),rho_ex))*_dt;
      rho_ex /= rho_ex.trace();
    }

    // Average state
    MatrixXcd rho = MatrixXcd::Zero(_dim_Hilbert_space, _dim_Hilbert_space);

    // Cycle on the ensemble members
    for (int i = 0; i < _N_states; ++i) {
      // Prints the trajectories
      if (verbose && i < _N_traj_print && _print_trajectory)
        traj << observable(projector(psi[i])) << " ";

      // Updates the average
      rho += projector(psi[i])/((double)_N_states);

      MatrixXcd R = J(projector(psi[i]),t) + 0.5*(C(projector(psi[i]), t)*projector(psi[i]) + projector(psi[i])*C(projector(psi[i]), t).adjoint());
      
      // Draws a random number and calculates whether the evolution is deterministic or via a jump
      double z = (double)rand()/((double)RAND_MAX);

      if (z < real(R.trace())*_dt) // Jump
        psi[i] = this->jump(R,z);
      else {// Free evolution
        MatrixXcd K = H(t) + 0.5*(C(projector(psi[i]), t).imag() - complex<double>(0.,1.)*(Gamma(t) + C(projector(psi[i]), t).real() ) );
        psi[i] -= K*psi[i]*complex<double>(0.,1.)*_dt;
      }
      psi[i] = psi[i].normalized();
    }
    // Storing the observable
    observables[n_observable] = observable(rho);
    n_observable++;

    if (verbose) {traj << endl;}
  }
  return observables;
}


// --- Running with all the copies
void roqj::run () {
  if (_verbose)
    print_info();

  ofstream params;
  params.open("params.txt");
  params << _N_copies << endl << _N_states << endl << _t_i << endl << _t_f << endl << _dt << endl << _print_trajectory << endl << _N_traj_print << endl << _dim_Hilbert_space;
  params.close();

  MatrixXd observables(_N_copies, _num_timesteps);
  for (int i = 0; i < _N_copies; ++i) {
    if (_verbose)
      cout << "Running copy " << i+1 << "/" << _N_copies << "...\n";
    observables.row(i) = run_single_iterations(i==0);
  }

  for (int i = 0; i < _num_timesteps; ++i) {
    _observable[i] = observables.col(i).mean();
    // Really there is no built-in method for the std deviation?
    _sigma_observable[i] = sqrt((observables.col(i).array() - _observable[i]).square().sum() / (observables.col(i).size() - 1));
  }
}

VectorXcd roqj::jump (const MatrixXcd &R, double z) const {
  ComplexEigenSolver<MatrixXcd> eigs;
  eigs.compute(R);
  VectorXcd eigval = eigs.eigenvalues();
  MatrixXcd eigvec = eigs.eigenvectors();
  //VectorXcd eigval;
  //MatrixXcd eigvec;

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
  return VectorXcd::Ones(_dim_Hilbert_space).normalized();
}

void roqj::print_info () const {
  cout << "\nRate Operator Quantum Jumps - running " << _N_copies << " copies.\n";
  cout << "\tNumber of states = " << _N_states << ", " << _dim_Hilbert_space << "-dimensional Hilbert space,\n";
  cout << "\tt_i = " << _t_i << ", t_f = " << _t_f << ", dt = " << _dt << ",\n";
  if (_print_trajectory)
    cout << "\tPrinting " << _N_traj_print << " trajectories.\n\n";
  else cout << endl;
}











// ------------------------- Qubit ROQJ class -------------------------
// --- Constructors
qubit_roqj::qubit_roqj (int N_states, double t_i, double t_f, double dt, int N_copies, bool print_trajectory, int N_traj_print, bool verbose) {
  srand(0);
  initialize(N_states, t_i, t_f, dt, N_copies, 2, print_trajectory, N_traj_print, verbose);
}

// -- Set initial state vector
void qubit_roqj::set_initial_state (const VectorXcd &psi) {
  if (psi.norm() == 0 || psi.size() != 2) {
    set_initial_state();
    return;
  }
  _initial_state = psi.normalized();
}

// Default initial state
void qubit_roqj::set_initial_state () {VectorXcd a; a << 1./sqrt(2.), 1./sqrt(2.); _initial_state = a;}

// --- Jump
VectorXcd qubit_roqj::jump (const MatrixXcd &R, double z) const {
  ComplexEigenSolver<MatrixXcd> eigs;
  eigs.compute(R);
  VectorXcd eigval = eigs.eigenvalues();
  MatrixXcd eigvec = eigs.eigenvectors();

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
  return VectorXcd::Ones(2).normalized();
}


// --- Run single iteration
VectorXd qubit_roqj::run_single_iterations (bool verbose) const {
  VectorXd observables(_num_timesteps);
  int n_observable = 0;

  // Allocating _N_states copies of the initial state
  std::vector<VectorXcd> psi;
  for (int i = 0; i <= _N_states; ++i)
    //psi[i] = _initial_state;
    psi.push_back(_initial_state);

  // Exact solution
  MatrixXcd rho_ex = MatrixXcd::Zero(2, 2);
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
      rho_ex /= rho_ex.trace();
    }

    // Average state
    MatrixXcd rho = MatrixXcd::Zero(2, 2);

    // Cycle on the ensemble members
    for (int i = 0; i < _N_states; ++i) {
      // Prints the trajectories
      if (verbose && i < _N_traj_print && _print_trajectory)
        traj << observable(projector(psi[i])) << " ";

      // Updates the average
      rho += projector(psi[i])/((double)_N_states);

      MatrixXcd R = J(projector(psi[i]),t) + 0.5*(C(projector(psi[i]), t)*projector(psi[i]) + projector(psi[i])*C(projector(psi[i]), t).transpose());
      
      // Draws a random number and calculates whether the evolution is deterministic or via a jump
      double z = (double)rand()/((double)RAND_MAX);

      if (z < real(R.trace())*_dt){ // Jump
        psi[i] = this->jump(R,z);
      }
      else {// Free evolution
        MatrixXcd K = H(t) + 0.5*(C(projector(psi[i]), t).imag() - complex<double>(0.,1.)*(Gamma(t) + C(projector(psi[i]), t).real() ) );
        psi[i] -= K*psi[i]*complex<double>(0.,1.)*_dt;
      }
      psi[i] = psi[i].normalized();
    }
    // Storing the observable
    observables[n_observable] = observable(rho);
    n_observable++;

    if (verbose) traj << endl;
  }
  return observables;
  return VectorXd::Zero(_num_timesteps);
}


// --- Running with all the copies
void qubit_roqj::run () {
  if (_verbose)
    print_info();

  ofstream params;
  params.open("params.txt");
  params << _N_copies << endl << _N_states << endl << _t_i << endl << _t_f << endl << _dt << endl << _print_trajectory << endl << _N_traj_print << endl << 2;
  params.close();

  MatrixXd matrix_observables = MatrixXd::Zero(_num_timesteps, _N_copies);
  for (int i = 0; i < _N_copies; ++i) {
    if (_verbose)
      cout << "Running copy " << i+1 << "/" << _N_copies << "...\n";
    VectorXd this_obs = run_single_iterations(i==0);
    for (int j = 0; j < _num_timesteps; ++j)
      matrix_observables(j,i) = this_obs(j);
  }

  for (int i = 0; i < _num_timesteps; ++i) {
    _observable[i] = matrix_observables.row(i).mean();
    _sigma_observable[i] = sqrt((matrix_observables.row(i).array() - _observable[i]).square().sum() / (matrix_observables.row(i).size() - 1));
  }
}











// ------------------------- ROQJ initially mixed -------------------------
roqj_mixed::roqj_mixed (int N_states, double t_i, double t_f, double dt, int N_copies, int dim_Hilbert_space, bool verbose, double threshold) {
  initialize(N_states, t_i, t_f, dt, N_copies, dim_Hilbert_space, false, 0, verbose, threshold);
}

// Setting the ensemble
void roqj_mixed::set_ensemble (const vector<pair<double, VectorXcd>> &ensemble) {
  _ensemble = ensemble;
  double sum_p = 0.;
  for (auto & elem : _ensemble) {
    if (elem.first < 0.)
      elem.first = -elem.first;
    sum_p += elem.first;
  }
  if (sum_p == 0.) {
    set_ensemble();
    return;
  }
  if (sum_p != 1.) {
    for (auto & elem : _ensemble)
      elem.first /= sum_p;
  }
}

void roqj_mixed::set_ensemble (const vector<double> &probabilities, const vector<VectorXcd> &states) {
  if (probabilities.size() != states.size()) {
    cerr << "Different number of probabilities and states -- EXIT\n";
    exit(EXIT_FAILURE);
  }
  _ensemble.clear();
  vector<pair<double, VectorXcd>> ens;
  for (int i = 0; i < probabilities.size(); ++i)
    ens.push_back(pair<double, VectorXcd>(probabilities[i], states[i]));
  set_ensemble(ens);
}

void roqj_mixed::set_ensemble () {
  VectorXcd init_state = VectorXcd::Ones(_dim_Hilbert_space).normalized();
  _ensemble.clear();
  _ensemble.push_back(pair<double, VectorXcd>(1., init_state));
}

void roqj_mixed::add_ensemble (const pair<double, VectorXcd> &state) {
  pair<double, VectorXcd> ens(state);
  if (ens.first == 0.) return;
  if (ens.first < 0.) ens.first = -ens.first;
  _ensemble.push_back(ens);

  double normalization = 0.;
  for (auto & i : _ensemble)
    normalization += i.first;
  if (normalization != 1.) {
    for (auto & i : _ensemble)
      i.first /= normalization;
  }
}

void roqj_mixed::add_ensemble (double prob, const VectorXcd &state) {
  add_ensemble(pair<double, VectorXcd>(prob, state));
}

// Runs the ROQJ for each state and takes the average state. Repeats it _N_copies time
void roqj_mixed::run () {
  if (_verbose) {
    print_info();
    print_ens();
  }

  ofstream params;
  params.open("params.txt");
  params << _N_copies << endl << _N_states << endl << _t_i << endl << _t_f << endl << _dt << endl << false << endl << 0 << endl << _dim_Hilbert_space;
  params.close();

  MatrixXd matrix_observables = MatrixXd::Zero(_num_timesteps, _N_copies);
  // Cycle on the copies
  for (int i = 0; i < _N_copies; ++i) {
    if (_verbose)
      cout << "Running copy " << i+1 << "/" << _N_copies << "...\n";
    
    // Cycle on the ensemble members'
    VectorXd obs = VectorXd::Zero(_num_timesteps);
    int n = 1;
    for (auto & elem : _ensemble) {
      set_initial_state(elem.second);
      if (_verbose)
        cout << "\tEnsemble member " << n++ << "/" << _ensemble.size() << "...\n";
      obs += elem.first * run_single_iterations(false);
      //obs = elem.first * run_single_iterations(false);
    }
    if (_verbose) cout << endl;

    for (int j = 0; j < _num_timesteps; ++j)
      matrix_observables(j,i) = obs(j);
  }

  for (int i = 0; i < _num_timesteps; ++i) {
    _observable[i] = matrix_observables.row(i).mean();
    _sigma_observable[i] = sqrt((matrix_observables.row(i).array() - _observable[i]).square().sum() / (matrix_observables.row(i).size() - 1));
  }
}


VectorXd roqj_mixed::get_exact_sol (string file_out) {
  ofstream out;
  bool verbose = file_out != "";
  if (verbose)
    out.open(file_out);

  MatrixXcd rho = MatrixXcd::Zero(_dim_Hilbert_space, _dim_Hilbert_space);
  for (auto & elem : _ensemble)
    rho += elem.first*projector(elem.second);

  VectorXd obs = VectorXd(_num_timesteps);
  int n = 0;

  for (double t = _t_i; t <= _t_f; t += _dt) {
    double this_obs = observable(rho);
    if (verbose)
      out << this_obs << endl;
    obs[n] = this_obs;
    n++;
    rho = rho + (-I*comm(H(t),rho) + J(rho,t) - 0.5*anticomm(Gamma(t),rho))*_dt;
  }

  return obs;
}

void roqj_mixed::print_ens () {
  cout << "Ensemble:\n";
  cout << "\tState\t| Probability\n";
  for (auto & elem : _ensemble) {
    cout << "----------------|----------------\n";
    cout << elem.second << "\t| " << elem.first << endl;
  }
  cout << endl;
  return;
}