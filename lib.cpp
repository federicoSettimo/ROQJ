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

cx_mat projector (const cx_vec &psi) {return psi*psi.t();}

// ------------------------- METHODS DEFINITIONS -------------------------
// --- Constructors
void roqj::initialize (int N_ensemble, double t_i, double t_f, double dt, int N_copies, int dim_Hilbert_space) {
  set_N_ensemble(N_ensemble);
  set_time(t_i, t_f, dt);
  set_N_copies(N_copies);
  set_dim_Hilbert_space(dim_Hilbert_space);
  set_initial_state();
}

roqj::roqj (int N_ensemble, double t_i, double t_f, double dt, int N_copies, int dim_Hilbert_space) {
  initialize(N_ensemble, t_i, t_f, dt, N_copies, dim_Hilbert_space);
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

void roqj::set_initial_state (const cx_vec &psi_i) {
  double n = arma::norm(psi_i);
  if (n == 0 || psi_i.size() != _dim_Hilbert_space) {
    set_initial_state();
    return;
  }
  _initial_state = psi_i;
  normalise(_initial_state);
}

void roqj::set_initial_state () {
  _initial_state = cx_vec(_dim_Hilbert_space, arma::fill::ones);
  normalise(_initial_state);
}

// --- Getters
int roqj::get_N_ensemble () const {return _N_ensemble;}
int roqj::get_N_copies () const {return _N_copies;}
int roqj::get_dim_Hilbert_space () const {return _dim_Hilbert_space;}
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
vec roqj::run_single_iterations () const {
  srand(time(NULL));

  vec observables(_num_timesteps);
  int n_observable = 0;

  // Allocating _N_ensemble copies of the initial state
  std::vector<cx_vec> psi(_N_ensemble);
  for (int i = 0; i <= _N_ensemble; ++i)
    psi[i] = _initial_state;

  // Time evolution
  for (double t = _t_i; t <= _t_f; t += _dt) {
    // Average state
    cx_mat rho(_dim_Hilbert_space, _dim_Hilbert_space, arma::fill::zeros);

    // Cycle on the ensemble members
    for (int i = 0; i < _N_ensemble; ++i) {
      // Updates the average
      rho += projector(psi[i])/((double)_N_ensemble);

      cx_mat R = J_prime(projector(psi[i]),t);
      
      // Draws a random number and calculates whether the evolution is deterministic or via a jump
      double z = (double)rand()/((double)RAND_MAX);

      if (z < real(arma::trace(R))) { // Jump
        cx_vec eigval;
        cx_mat eigvec;
        eig_gen(eigval, eigvec, R);

        // Chose in which eigenvalue perform the jump
        double sum_previous_eigs = 0.;
        bool already_jumped = false;
        for (int j = 0; j < _dim_Hilbert_space && !already_jumped; ++j) {
          // If z is in the j-th bin, it jumps to the j-th eigenstate
          if (z >= sum_previous_eigs && z < sum_previous_eigs + real(eigval[i])) {
            already_jumped = true;
            cx_vec post_jump_state(_dim_Hilbert_space);
            for (int k = 0; k < _dim_Hilbert_space; ++k) {
              // j-th eigenstate = j-th coloumn of the eigenvectors matrix
              post_jump_state[k] = eigvec(k,j);
            }
            psi[i] = post_jump_state;
          }
          sum_previous_eigs += real(eigval[i]);
        }
      }
      else { // Free evolution
        psi[i] -= K(rho, t)*psi[i]*complex<double>(0.,1.)*_dt;
        normalise(psi[i]);
      }
    }
    // Storing the observable
    observables[n_observable] = observable(rho);
    n_observable++;
  }
  return observables;
}


// --- Running with all the copies
void roqj::run () {
  cout << "Rate Operator Quantum Jumps - running " << _N_copies << " copies.\n";
  cout << "\tEnsemble size = " << _N_ensemble << ", " << _dim_Hilbert_space << "-dimensional Hilbert space\n";
  cout << "\tt_i = " << _t_i << ", t_f = " << _t_f << ", dt = " << _dt << endl << endl;

  ofstream params;
  params.open("params.txt");
  params << _N_copies << endl << _N_ensemble << endl << _t_i << endl << _t_f << endl << _dt << endl << _dim_Hilbert_space;

  for (int i = 0; i < _N_copies; ++i) {
    cout << "Running copy " << i+1 << "/" << _N_copies << "...\n";
    _observable += run_single_iterations()/((double)_N_copies);
  }
}