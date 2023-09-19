#include "roqj_state.h"
// ------------------------- FUNCTIONS DEFINITIONS -------------------------
bool isNormalized (const VectorXcd &psi) {return psi.norm() == 1;}

MatrixXcd comm (const MatrixXcd &A, const MatrixXcd &B) {return A*B-B*A;}

MatrixXcd anticomm (const MatrixXcd &A, const MatrixXcd &B) {return A*B+B*A;}

MatrixXcd projector (const VectorXcd &psi) {return psi*psi.adjoint();}

MatrixXcd BlochToMatrix (double x, double y, double z) {
  double r = sqrt(x*x + y*y + z*z);
  if (r > 1.) {x /= r; y /= r; z /= r;}
  return .5*(id + x*sigma_x + y*sigma_y + z*sigma_z);
}

MatrixXcd tr_1(const MatrixXcd &rho) {
	MatrixXcd A(2,2);
  A << rho(0,0) + rho(2,2), rho(0,1) + rho(2,3), rho(1,0) + rho(3,2), rho(1,1) + rho(3,3);
	return A;
}

MatrixXcd tr_2(const MatrixXcd &rho) {
	MatrixXcd A(2,2);
  A << rho(0,0) + rho(1,1), rho(0,2) + rho(1,3), rho(2,0) + rho(3,1), rho(2,2) + rho(3,3);
  return A;
}

MatrixXcd tens (const MatrixXcd &A, const MatrixXcd &B) {
  MatrixXcd C;
  C = MatrixXcd::Zero(4,4);
  C.topLeftCorner(2,2) = A(0,0)*B;
  C.topRightCorner(2,2) = A(0,1)*B;
  C.bottomLeftCorner(2,2) = A(1,0)*B;
  C.bottomRightCorner(2,2) = A(1,1)*B;
  return C;
}


VectorXcd tens_state (const Vector2cd &psi1, const Vector2cd &psi2) {
  VectorXcd psi;
  psi = VectorXcd::Zero(4);
  psi(0) = psi1(0)*psi2(0);
  psi(1) = psi1(0)*psi2(1);
  psi(2) = psi1(1)*psi2(0);
  psi(3) = psi1(1)*psi2(1);
  return psi;
}

double entropy (const Matrix2cd &rho) {
  ComplexEigenSolver<Matrix2cd> eigs;
  eigs.compute(rho);
  double p = real(eigs.eigenvalues()(0));
  return -p*log2(p) - (1.-p)*log2(1.-p);
}




// ------------------------- METHODS DEFINITIONS -------------------------
// ------------------------- ROQJ class -------------------------
// --- Constructors
void roqj::initialize (int N_states, double t_i, double t_f, double dt, int N_copies, int dim_Hilbert_space, bool print_trajectory, int N_traj_print, bool verbose, double threshold) {
  set_N_states(N_states);
  set_time(t_i, t_f, dt);
  set_N_copies(N_copies);
  set_dim_Hilbert_space(dim_Hilbert_space);
  set_print_traj(print_trajectory, N_traj_print);
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
  _num_timesteps = (int)((_t_f - _t_i)/_dt) + 1;
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
 set_print_traj(_print_trajectory, N_traj_print);
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

void roqj::set_print_traj (bool print, int N_traj_print) {
  if (print && N_traj_print > 0) {
    _print_trajectory = true;
    _N_traj_print = N_traj_print;
  }
  else {
    _print_trajectory = false;
    _N_traj_print = 0;
  }
}

void roqj::set_verbose (bool verbose) {_verbose = verbose;}

void roqj::set_threshold (double threshold) {
  if (threshold < 0.)
    threshold = -threshold;
  if (threshold > 1.)
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

VectorXd roqj::get_observable (string file_out) const {
  if (file_out == "")
    return _observable;
  ofstream out;
  out.open(file_out);
  for (int i = 0; i < _num_timesteps; ++i)
    out << _observable[i] << endl;
  return _observable;
}

VectorXd roqj::get_error_observable (string file_out) const {
  if (file_out == "")
    return _sigma_observable;
  ofstream out;
  out.open(file_out);
  for (int i = 0; i < _num_timesteps; ++i)
    out << _sigma_observable[i] << endl;
  return _sigma_observable;
}


VectorXd roqj::get_det_trajectory (string file_out) const {
  VectorXcd psi = _initial_state;
  VectorXd traj(_num_timesteps);
  int i = 0;
  ofstream out;
  if (file_out != "")
    out.open(file_out);
  for (double t = _t_i; t <= _t_f; t += _dt) {
    MatrixXcd K = H(t) - 0.5*I*Gamma(t);
    psi -=  K*psi*I*_dt + .5*_dt*Phi(psi,t,false);
    psi = psi.normalized();
    double x = observable(projector(psi));
    traj[i] = x;
    ++i;
    if (file_out != "")
      out << x << endl;
  }
  return traj;
}


// --- Run single iteration
VectorXd roqj::run_single_iterations (bool verbose) const {
  VectorXd observables(_num_timesteps);
  int n_observable = 0;

  // Allocating _N_states copies of the initial state
  std::vector<VectorXcd> psi;
  vector<bool> jumped(_N_states);
  for (int i = 0; i <= _N_states; ++i) {
    psi.push_back(_initial_state);
    jumped[i] = false;
  }

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
      if (real(psi[i](0)) < 0. || (real(psi[i](0)) < .001 && real(psi[i](1)) < 0.)) psi[i] *= -1.;
      // Prints the trajectories
      if (verbose && i < _N_traj_print && _print_trajectory)
        traj << observable(projector(psi[i])) << " ";

      // Updates the average
      rho += projector(psi[i])/((double)_N_states);

      VectorXcd Phi_i = Phi(psi[i],t,jumped[i]);
      MatrixXcd R = J(projector(psi[i]),t) + 0.5*(Phi_i*(psi[i].adjoint()) + psi[i]*(Phi_i.adjoint()));
      
      // Draws a random number and calculates whether the evolution is deterministic or via a jump
      double z = (double)rand()/((double)RAND_MAX);

      if (z < real(R.trace())*_dt || real(R.trace()) < 0) {// Jump
        psi[i] = this->jump(R,z,psi[i]);
        jumped[i] = true;
      }
      else {// Free evolution
        MatrixXcd K = H(t) - 0.5*I*Gamma(t);
        psi[i] = psi[i] - K*psi[i]*I*_dt - .5*_dt*Phi_i;
      }
      psi[i] = psi[i].normalized();
      for (int j = 0; j < _dim_Hilbert_space; ++j) {
        if (abs(psi[i](j)) < _threshold)
          psi[i](j) = 0.;
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

VectorXcd roqj::jump (const MatrixXcd &R, double z, const VectorXcd &psi) const {
  ComplexEigenSolver<MatrixXcd> eigs;
  eigs.compute(R);
  VectorXcd eigval = eigs.eigenvalues();
  MatrixXcd eigvec = eigs.eigenvectors();

  // Chose in which eigenvalue perform the jump
  double sum_previous_eigs = 0.;
  for (int j = 0; j < _dim_Hilbert_space; ++j) {
    if (real(eigval[j]) < -_threshold) {
      cerr << "Negative rate - reverse jump. NOT IMPLEMENTED\n";
      cout << "State: (";
      for (int i = 0; i < _dim_Hilbert_space; ++i) {
        double re = real(psi(i)), im = imag(psi(i));
        if (re == 0 && im == 0)
          cout << 0;
        else if (re == 0)
          cout << im << " I";
        else if (im == 0)
          cout << re;
        else {
          cout << re;
          if (im > 0)
            cout << " + ";
          else cout << " - ";
          cout << abs(im) << " I";
        }

        if (i == _dim_Hilbert_space - 1)
          cout << ")";
        else cout << ",  ";
      }
      cout << "\nEigenvalues: ";
      for (int i = 0; i < _dim_Hilbert_space; ++i) 
        cout << real(eigval[i]) << ", ";
      cout << endl;
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
qubit_roqj::qubit_roqj (int N_states, double t_i, double t_f, double dt, int N_copies, bool print_trajectory, int N_traj_print, bool verbose, double threshold) {
  srand(0);
  initialize(N_states, t_i, t_f, dt, N_copies, 2, print_trajectory, N_traj_print, verbose, threshold);
}

// -- Set initial state vector
void qubit_roqj::set_initial_state (const VectorXcd &psi) {
  if (psi.norm() == 0 || psi.size() != 2) {
    set_initial_state();
    return;
  }
  _initial_state = psi.normalized();
  return;
}

// Default initial state
void qubit_roqj::set_initial_state () {VectorXcd a; a << 1./sqrt(2.), 1./sqrt(2.); _initial_state = a;}

// --- Jump
VectorXcd qubit_roqj::jump (const MatrixXcd &R, double z, const VectorXcd &psi) const {
  ComplexEigenSolver<MatrixXcd> eigs;
  eigs.compute(R);
  VectorXcd eigval = eigs.eigenvalues();
  MatrixXcd eigvec = eigs.eigenvectors();

  VectorXcd phi1 = eigvec.col(0), phi2 = eigvec.col(1);
  if (real(phi1(0)) < 0. || (real(phi1(0)) < .0001 && real(phi1(1)) < 0.)) phi1 = -phi1;
  if (real(phi2(0)) < 0. || (real(phi2(0)) < .0001 && real(phi2(1)) < 0.)) phi2 = -phi2;

  double lambda1 = real(eigval[0]), lambda2 = real(eigval[1]), pjump1 = lambda1*_dt;
  if (lambda1 >= -_threshold && lambda2 >= -_threshold) {// Normal jump
    // With probability pjump1, it jumps to the first eigenstate of R
    if (z <= pjump1)
      return phi1;
    else return phi2;
  }
  else {// Reverse jump ----- Not implemented??
    cerr << "Negative rate - reverse jump. NOT IMPLEMENTED\n";
    cout << "Eigenvalues: " << lambda1 << ", " << lambda2 << endl;
    cout << "State: " << psi.transpose() << endl;
    exit(EXIT_FAILURE);
  }
  return VectorXcd::Ones(2).normalized();
}