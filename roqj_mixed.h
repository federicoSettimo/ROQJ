/*
  ROQJ class for initially mixed states
*/
#ifndef _ROQJ_MIXED_H_
#define _ROQJ_MIXED_H_

#include "roqj.h"
#include <utility>

// ------------------------- ROQJ for mixed initial state -------------------------
class roqj_mixed : public roqj{
protected:
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

  // Printing the ensemble
  void print_ens ();

  // Ad a state to the ensemble
  void add_ensemble (const pair<double, VectorXcd> &state);
  void add_ensemble (double prob, const VectorXcd &state);

  // Runs the ROQJ for each state and takes the average state. Repeats it _N_copies time
  void run ();

  // Returns the exact value for the observable
  VectorXd get_exact_sol (string file_out = "");
};

#endif