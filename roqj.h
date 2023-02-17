// For a reasonbale execution time, use the flags
// -ffast-math -fno-math-errno -O3
// both for the library and for the program

/*
  Class producing the ROQJ trajectories.
  Required an external definition of the operators defining the dynamics (H, J, Gamma), of C defining the diferent ROs, and of the chosen observble
*/
#ifndef _ROQJ_H_
#define _ROQJ_H_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>
#include <vector>
#include <utility>
#include <cstring>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace Eigen;

// Default values
const int N_states_default = 10000, N_copies_default = 1, dim_Hilbert_space_default = 2, N_traj_print_default = 5;
const double t_i_default = 0., t_f_default = 10., dt_default = 0., threshold_default = 1.e-10;

// External functions needed: H(t), J(rho, t), Gamma(t), C(t), observable(rho)
extern MatrixXcd H (double t);
extern MatrixXcd J (const MatrixXcd &rho, double t);
extern MatrixXcd Gamma (double t);
extern MatrixXcd C (const MatrixXcd &rho, double t);
extern double observable (const MatrixXcd &rho);

// ------------------------- ROQJ class -------------------------
class roqj {
protected:
  int _N_states, _N_copies, _dim_Hilbert_space, _num_timesteps, _N_traj_print;
  bool _print_trajectory, _verbose;
  double _t_i, _t_f, _dt, _threshold;
  VectorXd _observable, _sigma_observable;
  VectorXcd _initial_state;
public:
  /* 
    Parameters: (int) ensemble size, (double) intial time, (double) final time, (double) dt, (int) number of copies, (int) dim Hilbert space, (bool) print trajectory, (int) number of trajectories to print, (bool) verbose, (double) threshold for negativity
    Default values: N_states = 10000, t_i = 0, t_f = 10, dt = t_f/10000, N_copies = 1, dim_Hilbert_space = 2, print_trajectory = true, N_traj_print = 3, verbose = true, threshold = 1e-20
  */
  roqj (int N_states = N_states_default, double t_i = t_i_default, double t_f = t_f_default, double dt = dt_default, int N_copies = N_copies_default, int dim_Hilbert_space = dim_Hilbert_space_default, bool print_trajectory = true, int N_traj_print = N_traj_print_default, bool verbose = true, double threshold = threshold_default);
  void initialize (int N_states = N_states_default, double t_i = t_i_default, double t_f = t_f_default, double dt = dt_default, int N_copies = N_copies_default, int dim_Hilbert_space = dim_Hilbert_space_default, bool print_trajectory = true, int N_traj_print = N_traj_print_default, bool verbose = true, double threshold = threshold_default);


  // Setting the initial state. If psi_i is not a dim_Hilbert_space-dimensional VectorXdtor, default initializer
  void set_initial_state (const VectorXcd &psi_i);
  // Setting the initial state. Default is (|1>+...+|dim_Hilbert_space>)/sqrt(dim_Hilbert_space). Returns false if psi_i is not a dim_Hilbert_space-dimensional VectorXdtor
  void set_initial_state ();


  // runs N_copies times the single iteration. Updates the observable and its error
  void run ();


  // One single iteration; to be run N_copies times. verbose = true: prints exact sol and trajectories
  VectorXd run_single_iterations (bool verbose = true) const;


  // Performs the jump process
  VectorXcd jump (const MatrixXcd &R, double z) const;


  // Displays the info on the runs
  void print_info () const;


  // Setters
  void set_N_states (int N_states = N_states_default);
  void set_N_copies (int N_states = N_copies_default);
  void set_t_i (double t_i = t_i_default);
  void set_t_f (double t_f = t_f_default);
  void set_dt (double dt = dt_default);
  void set_time (double t_i = t_i_default, double t_f = t_f_default, double dt = dt_default);
  void set_dim_Hilbert_space (int dim_Hilbert_space = dim_Hilbert_space_default);
  void set_N_traj_print (int N_traj_print = N_traj_print_default);
  void set_print_traj (bool print = true);
  void set_verbose (bool verbose = true);
  void set_threshold (double threshold = threshold_default);


  // Getters
  int get_N_states () const;
  int get_N_copies () const;
  int get_dim_Hilbert_space () const;
  int get_N_traj_print () const;
  double get_t_i () const;
  double get_t_f () const;
  double get_dt () const;
  double get_threshold () const;
  VectorXcd get_initial_state () const;
  
  // Returns the values of the observable
  VectorXd get_observable () const;
  // Prints the values of the observable in file_out
  VectorXd get_observable (string file_out) const;

  // Returns the errors of the observable
  VectorXd get_error_observable () const;
  // Prints the errors of the observable in file_out
  VectorXd get_error_observable (string file_out) const;
};





// ------------------------- Qubit ROQJ class -------------------------
class qubit_roqj:public roqj {
protected:
  Vector2cd _initial_state;
public:
  /* 
    Parameters: (int) ensemble size, (double) intial time, (double) final time, (double) dt, (int) number of copies, (bool) print trajectory, (int) number of trajectories to print, (bool) verbose
    Default values: N_states = 10000, t_i = 0, t_f = 10, dt = t_f/10000, N_copies = 1, print_trajectory = true, N_traj_print = 3, verbose = true
  */
  qubit_roqj (int N_states = N_states_default, double t_i = t_i_default, double t_f = t_f_default, double dt = dt_default, int N_copies = N_copies_default, bool print_trajectory = true, int N_traj_print = N_traj_print_default, bool verbose = true);

  // Setting initial state with a 2-d VectorXdtor
  void set_initial_state (const VectorXcd &psi);
  // Default initializer - Id/2
  void set_initial_state ();

  // Performs the jump with only 2 possible channels
  VectorXcd jump (const MatrixXcd &R, double z) const;

  // Runs with the 2-channel jumps
  void run ();

  // Single iteration with the 2-channel jumps
  VectorXd run_single_iterations (bool verbose = true) const;
};





// ------------------------- Qubit ROQJ with 2 fixed post-jump states -------------------------
class roqj_mixed {
protected:
  // Parameters of the roqj
  int _N_states, _N_copies, _dim_Hilbert_space, _num_timesteps;
  bool _verbose;
  double _t_i, _t_f, _dt, _threshold;
  VectorXd _observable, _sigma_observable;

  // ROQJ
  roqj _solver;

  // Ensemble
  vector<pair<double, VectorXcd>> _ensemble;

public:
  /* 
    Parameters: (int) ensemble size, (double) intial time, (double) final time, (double) dt, (int) number of copies, (int) dim Hilbert space, (bool) verbose, (double) threshold for negativity
    Default values: N_states = 10000, t_i = 0, t_f = 10, dt = t_f/10000, N_copies = 1, dim_Hilbert_space = 2, verbose = true, threshold = 1e-20
  */
  roqj_mixed (int N_states = N_states_default, double t_i = t_i_default, double t_f = t_f_default, double dt = dt_default, int N_copies = N_copies_default, int dim_Hilbert_space = dim_Hilbert_space_default, bool verbose = true, double threshold = threshold_default);

  // Setting the ensemble
  void set_ensemble (const vector<pair<double, VectorXcd>> &ensemble);
  void set_ensemble (const vector<double> &probabilities, const vector<VectorXcd> &states);
  void set_ensemble ();

  // Runs the ROQJ for each state and takes the average state. Repeats it _N_copies time
  void run ();

  // Runs one single iteration for each member of the ensebles, returns the observable for the average at each time
  VectorXd run_single_iteration ();

  // Returns the values of the observable
  VectorXd get_observable () const;
  // Prints the values of the observable in file_out
  VectorXd get_observable (string file_out) const;

  // Returns the errors of the observable
  VectorXd get_error_observable () const;
  // Prints the errors of the observable in file_out
  VectorXd get_error_observable (string file_out) const;

  // Returns the exact value for the observable
  VectorXd get_exact_sol (string file_out = "");
};





// ------------------------- FUNCTIONS -------------------------
// Checks whether the VectorXdtor is normalized
bool isNormalized (const VectorXcd &psi);

// Commutator
MatrixXcd comm (const MatrixXcd &A, const MatrixXcd &B);

// Anticommutator
MatrixXcd anticomm (const MatrixXcd &A, const MatrixXcd &B);

// Projector |psi><psi|
MatrixXcd projector (const VectorXcd &psi);

// Density operator from its Bloch VectorXdtor representation
MatrixXcd BlochToMatrix (double x, double y, double z);





// ------------------------- PAULI MATRICES -------------------------
static complex<double> I(0,1), one(1,0);
static Eigen::MatrixXcd sigma_x {{0,1},{1,0}};
static Eigen::MatrixXcd sigma_y {{0,-I},{I,0}};
static Eigen::MatrixXcd sigma_z {{1,0},{0,-1}};
static Eigen::MatrixXcd sigma_p {{0,1},{0,0}};
static Eigen::MatrixXcd sigma_m {{0,0},{1,0}};
static Eigen::MatrixXcd id {{1,0},{0,1}};

#endif