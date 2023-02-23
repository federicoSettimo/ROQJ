/*
  ROQJ class for the special case in which the effective ensemble
  is finite dimensional.
  The dynamics is calculated by keeping track of the populations of the
  members of the effective ensemble.
*/
#ifndef _ROQJ_POP_H_
#define _ROQJ_POP_H_

#include "roqj.h"

extern const int N_states_default, N_copies_default, dim_Hilbert_space_default, N_traj_print_default;
extern const double t_i_default, t_f_default, dt_default, threshold_default;

// ------------------------- Qubit ROQJ pop class -------------------------
class qubit_roqj_pop : public qubit_roqj {
protected:
  Vector2cd _eig_1, _eig_2;
public:
  /* 
    Parameters: (int) ensemble size, (double) intial time, (double) final time, (double) dt, (int) number of copies, (bool) verbose
    Default values: N_states = 10000, t_i = 0, t_f = 10, dt = t_f/10000, N_copies = 1, print_trajectory = true, N_traj_print = 3, verbose = true
  */
  qubit_roqj_pop (int N_states = N_states_default, double t_i = t_i_default, double t_f = t_f_default, double dt = dt_default, int N_copies = N_copies_default, bool verbose = true, double threshold = threshold_default);

  // Runs with the 2-channel jumps
  void run ();

  // Single iteration with the 2-channel jumps
  VectorXd run_single_iterations (bool verbose = true) const;
};


// ------------------------- Qubit ROQJ pop initially mixed class -------------------------
class qubit_roqj_pop_mixed : public qubit_roqj {
protected:
  vector<pair<double, Vector2cd>> _ensemble;
  qubit_roqj_pop _roqj;
public:
  /* 
    Parameters: (int) ensemble size, (double) intial time, (double) final time, (double) dt, (int) number of copies, (bool) verbose
    Default values: N_states = 10000, t_i = 0, t_f = 10, dt = t_f/10000, N_copies = 1, print_trajectory = true, N_traj_print = 3, verbose = true
  */
  qubit_roqj_pop_mixed (int N_states = N_states_default, double t_i = t_i_default, double t_f = t_f_default, double dt = dt_default, int N_copies = N_copies_default, bool verbose = true, double threshold = threshold_default);

  void set_ensemble (const vector<pair<double, Vector2cd>> &ensemble);
  void set_ensemble (const vector<double> &probabilities, const vector<Vector2cd> &states);
  void set_ensemble ();

  // Printing the ensemble
  void print_ens ();

  // Ad a state to the ensemble
  void add_ensemble (const pair<double, Vector2cd> &state);
  void add_ensemble (double prob, const Vector2cd &state);

  // Runs the ROQJ for each state and takes the average state. Repeats it _N_copies time
  void run ();

  // Returns the exact value for the observable
  VectorXd get_exact_sol (string file_out = "");
};




#endif