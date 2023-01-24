/*
  Class producing the ROQJ trajectories.
  Required an external definition of the operators defining the dynamic (H, J, Gamma), of C defining the diferent ROs, and of the chosen observble
*/
#ifndef _ROQJ_H_
#define _ROQJ_H_

#include <armadillo>
#include <complex>
#include <vector>
#include <cstring>
#include <fstream>
#include <ctime>

/*
  TO DO:
  - add the possibility to also print some of the trajectories
    In run_single_iteration add (int)n_iteartion and (bool)print
*/

using namespace std;
using namespace arma;

// Default values
const int N_ensemble_default = 10000, N_copies_default = 1, dim_Hilbert_space_default = 2, N_traj_print_default = 5;
const double t_i_default = 0., t_f_default = 10., dt_default = 0.;

// External functions needed: H(t), J(rho, t), Gamma(t), C(t), observable(rho)
extern cx_mat H (double t);
extern cx_mat J (const cx_mat &rho, double t);
extern cx_mat Gamma (double t);
extern cx_mat C (const cx_mat &rho, double t);
extern double observable (const cx_mat &rho);

// ------------------------- ROQJ class -------------------------
class roqj {
private:
  int _N_ensemble, _N_copies, _dim_Hilbert_space, _num_timesteps, _N_traj_print;
  bool _print_trajectory;
  double _t_i, _t_f, _dt;
  vec _observable, _sigma_observable;
  cx_vec _initial_state;
public:
  /* 
    Parameters: (int) ensemble size, (double) intial time, (double) final time, (double) dt, (int) number of copies, (int) dim Hilbert space, (bool) print trajectory, (int) number of trajectories to print
    Default values: N_ensemble = 10000, t_i = 0, t_f = 10, dt = t_f/10000, N_copies = 1, dim_Hilbert_space = 2, print_trajectory = true, N_traj_print = 3
  */
  roqj (int N_ensemble = N_ensemble_default, double t_i = t_i_default, double t_f = t_f_default, double dt = dt_default, int N_copies = N_copies_default, int dim_Hilbert_space = dim_Hilbert_space_default, bool print_trajectory = true, int N_traj_print = N_traj_print_default);
  void initialize (int N_ensemble = N_ensemble_default, double t_i = t_i_default, double t_f = t_f_default, double dt = dt_default, int N_copies = N_copies_default, int dim_Hilbert_space = dim_Hilbert_space_default, bool print_trajectory = true, int N_traj_print = N_traj_print_default);


  // Setting the initial state.If psi_i is not a dim_Hilbert_space-dimensional vector, default initializer
  void set_initial_state (const cx_vec &psi_i);
  // Setting the initial state. Default is (|1>+...+|dim_Hilbert_space>)/sqrt(dim_Hilbert_space). Returns false if psi_i is not a dim_Hilbert_space-dimensional vector
  void set_initial_state ();


  // runs N_copies times the single iteration. Updates the observable and its error
  void run ();


  // One single iteration; to be run N_copies times. verbose = true: prints exact sol and trajectories
  vec run_single_iterations (bool verbose = true) const;

  void reset ();

  // Setters
  void set_N_ensemble (int N_ensemble = N_ensemble_default);
  void set_N_copies (int N_ensemble = N_copies_default);
  void set_t_i (double t_i = t_i_default);
  void set_t_f (double t_f = t_f_default);
  void set_dt (double dt = dt_default);
  void set_time (double t_i = t_i_default, double t_f = t_f_default, double dt = dt_default);
  void set_dim_Hilbert_space (int dim_Hilbert_space = dim_Hilbert_space_default);
  void set_N_traj_print (int N_traj_print = N_traj_print_default);


  // Getters
  int get_N_ensemble () const;
  int get_N_copies () const;
  int get_dim_Hilbert_space () const;
  int get_N_traj_print () const;
  double get_t_i () const;
  double get_t_f () const;
  double get_dt () const;
  cx_vec get_initial_state () const;
  
  // Returns the values of the observable
  vec get_observable () const;
  // Prints the values of the observable in file_out
  vec get_observable (string file_out) const;

  // Returns the errors of the observable
  vec get_error_observable () const;
  // Prints the errors of the observable in file_out
  vec get_error_observable (string file_out) const;
};

// ------------------------- FUNCTIONS -------------------------
// Checks whether the vector is normalised
bool isNormalised (const cx_vec &psi);

// Commutator
cx_mat comm (const cx_mat &A, const cx_mat &B);

// Anticommutator
cx_mat anticomm (const cx_mat &A, const cx_mat &B);

// Projector |psi><psi|
cx_mat projector (const cx_vec &psi);

#endif