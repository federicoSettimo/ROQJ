#ifndef _ROQJ_H_
#define _ROQJ_H_

#include <armadillo>
#include <complex>
#include <vector>
#include <cstring>
#include <fstream>
#include <ctime>

const int N_ensemble_default = 10000, N_copies_default = 1, dim_Hilbert_space_default = 2;
const double t_i_default = 0., t_f_default = 10., dt_default = 0.;

// External functions needed: H(t), J(rho, t), Gamma(t), C(t), observable(rho)
extern cx_mat H (double t);
extern cx_mat J (const cx_mat &rho, double t);
extern cx_mat Gamma (double t);
extern cx_mat C (double t);
extern double observable (const cx_mat &rho);

class roqj {
private:
  int _N_ensemble, _N_copies, _dim_Hilbert_space, _num_observable;
  double _t_i, _t_f, _dt, _dt_observable;
  vec _observable, _sigma_observable; // Observable and its error
  cx_vec _initial_state;
public:
  /* 
    Parameters: (int) ensemble size, (double) intial time, (double) final time, (double) dt, (int) number of copies, (double) dt for calculating the observables
    Default values: N_ensemble = 10000, t_i = 0, t_f = 10, dt = t_f/10000, N_copies = 1
  */
  roqj (int N_ensemble = N_ensemble_default, double t_i = t_i_default, double t_f = t_f_default, double dt = dt_default, int N_copies = N_copies_default, int dim_Hilbert_space = dim_Hilbert_space_default, double dt_observable = dt_default);
  void initialize (int N_ensemble = N_ensemble_default, double t_i = t_i_default, double t_f = t_f_default, double dt = dt_default, int N_copies = N_copies_default, int dim_Hilbert_space = dim_Hilbert_space_default, double dt_observable = dt_default);


  // Setting the initial state.If psi_i is not a dim_Hilbert_space-dimensional vector, default initializer
  void set_initial_state (const cx_vec &psi_i);
  // Setting the initial state. Default is (|1>+...+|dim_Hilbert_space>)/sqrt(dim_Hilbert_space). Returns false if psi_i is not a dim_Hilbert_space-dimensional vector
  void set_initial_state ();


  // runs N_copies times the single iteration. Updates the observable and its error
  void run ();


  // One single iteration; to be run N_copies times
  vec run_single_iterations () const;


  // Setters
  void set_N_ensemble (int N_ensemble = N_ensemble_default);
  void set_N_copies (int N_ensemble = N_copies_default);
  void set_t_i (double t_i = t_i_default);
  void set_t_f (double t_f = t_f_default);
  void set_dt (double dt = dt_default);
  void set_time (double t_i = t_i_default, double t_f = t_f_default, double dt = dt_default);
  void set_dim_Hilbert_space (int dim_Hilbert_space = dim_Hilbert_space_default);
  void set_dt_observable (double dt_observable = dt_default);


  // Getters
  int get_N_ensemble () const;
  int get_N_copies () const;
  int get_dim_Hilbert_space () const;
  double get_t_i () const;
  double get_t_f () const;
  double get_dt () const;
  double get_dt_observable () const;
  
  // Returns the values of the observable
  vec get_observable () const;
  // Prints the values of the observable in file_out
  vec get_observable (string file_out) const;

  // Returns the errors of the observable
  vec get_error_observable () const;
  // Prints the errors of the observable in file_out
  vec get_error_observable (string file_out) const;
};

// Normalizes the state vector
void normalize (cx_vec &psi);

// Checks whether the vector is normalized
bool isNormalized (const cx_vec &psi);

// Commutator
cx_mat comm (const cx_mat &A, const cx_mat &B);

// ------------------------- FUNCTIONS DEFINITIONS -------------------------
void normalize (cx_vec &psi) {
  double n = psi.norm();
  if (n > 0.) psi = psi/n;
}

bool isNormalized (const cx_vec &psi) {return psi.norm() == 1;}

cx_mat comm (const cx_mat &A, const cx_mat &B) {return A*B-B*A;}

// ------------------------- METHODS DEFINITIONS -------------------------
// Constructors
void roqj::initialize (int N_ensemble, double t_i, double t_f, double dt, int N_copies, int dim_Hilbert_space, double dt_observable) {
  set_N_ensemble(N_ensemble);
  set_time(t_i, t_f, dt);
  set_N_copies(N_copies);
  set_dim_Hilbert_space(dim_Hilbert_space);
  set_dt_observable(dt_observable);
}

roqj::roqj (int N_ensemble, double t_i, double t_f, double dt, int N_copies, int dim_Hilbert_space, double dt_observable) {
  initialize(N_ensemble, t_i, t_f, dt, N_copies, dim_Hilbert_space, dt_observable);
}

// Setter
void roqj::set_N_ensemble (int N_ensemble = N_ensemble_default) {
  if (N_ensemble <= 0) _N_ensemble = N_ensemble_default;
  else _N_ensemble = N_ensemble;
}

void roqj::set_N_copies (int N_copies = N_copies_default) {
  if (N_copies <= 0) _N_copies = N_copies_default;
  else _N_copies = N_copies;
}

void roqj::set_t_i (double t_i = t_i_default) {
  if (t_i >= _t_f || t_i < 0.) _t_i = t_f_default;
  else _t_i = t_i;
}

void roqj::set_t_f (double t_f = t_f_default) {
  if (t_f <= _t_i) {
    if (t_f_default > _t_i) _t_f = t_f_default;
    else _t_f = 10.*_t_i;
  }
  else _t_f = t_f;
}

void roqj::set_dt (double dt = dt_default) {
  if (dt == dt_default || dt < 0. || dt >= _t_f - _t_i) _dt = dt = _t_f/10000.;
  else _dt = dt;
}
void roqj::set_time (double t_i = t_i_default, double t_f = t_f_default, double dt = dt_default) {
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

void roqj::set_dim_Hilbert_space (int dim_Hilbert_space = dim_Hilbert_space_default) {
  if (dim_Hilbert_space < 2) _dim_Hilbert_space = dim_Hilbert_space_default;
  else _dim_Hilbert_space = dim_Hilbert_space;
}

void roqj::set_dt_observable (double dt_observable = dt_default) {
  if (dt_observable <= _dt) _dt_observable = dt_default;
  else _dt_observable = dt_observable;
  _num_observable = (int)(_t_f - _t_i)/_dt_observable;
}

void roqj::set_initial_state (const cx_vec &psi_i) {
  n = psi_i.norm();
  if (n == 0 || size(psi_i) != _dim_Hilbert_space) {
    set_initial_state();
    return;
  }
  _initial_state = psi_i;
  normalize(_initial_state);
}

void roqj::set_initial_state () {
  _initial_state = cx_vec(_dim_Hilbert_space, fill::ones);
  normalize(_initial_state);
}

// Setters
int roqj::get_N_ensemble () const {return _N_ensemble;}
int roqj::get_N_copies () const {return _N_copies;}
int roqj::get_dim_Hilbert_space () const {return _dim_Hilbert_space;}
double roqj::get_t_i () const {return _t_i;}
double roqj::get_t_f () const {return _t_f;}
double roqj::get_dt () const {return _dt;}
double roqj::get_dt_observable () const {return _dt_observable;}
vec roqj::get_observable () const {return _observable;}
vec roqj::get_error_observable () const {return _sigma_observable;}

vec roqj::get_observable (string file_out) const {
  ofstream out;
  out.open(file_out);
  for (int i = 0; i < _num_observable; ++i)
    out << _observable[i] << endl;
}

vec roqj::get_error_observable (string file_out) const {
  ofstream out;
  out.open(file_out);
  for (int i = 0; i < _num_observable; ++i)
    out << _error_observable[i] << endl;
}


// Run single iteration
vec roqj::run_single_iterations () const {
  srand(time(NULL));

  vector<vec> psi(_N_ensemble, _initial_state);
}

// allocates a vector<vec> with size N_ensemble
// Returns the value of the observable

#endif _ROQJ_H_