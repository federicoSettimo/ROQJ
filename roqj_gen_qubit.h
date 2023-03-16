/*
  Class using the generic qubit ROQJ method. The main difference from the roqj class
  is that before the eigs of R were calculted only when performing the jumps.
  Now, since C already requires to calculate the eigs, they are calculated there and using the same eigenstates
  with modified eigenvalues for the jump method
*/

#ifndef _ROQJ_GEN_QUBIT_H_
#define _ROQJ_GEN_QUBIT_H_

#include "roqj.h"

class gen_qubit_roqj : public qubit_roqj {
public:
  /* 
    Parameters: (int) ensemble size, (double) intial time, (double) final time, (double) dt, (int) number of copies, (bool) print trajectory, (int) number of trajectories to print, (bool) verbose
    Default values: N_states = 10000, t_i = 0, t_f = 10, dt = t_f/10000, N_copies = 1, print_trajectory = true, N_traj_print = 3, verbose = true
  */
  gen_qubit_roqj (int N_states = N_states_default, double t_i = t_i_default, double t_f = t_f_default, double dt = dt_default, int N_copies = N_copies_default, bool print_trajectory = true, int N_traj_print = N_traj_print_default, bool verbose = true, double threshold = threshold_default);

  // Performs the jump with only 2 possible channels
  Vector2cd jump (const Vector2cd &phi1, const Vector2cd &phi2, double lambda1, double lambda2, double z)  const;

  // Runs with the 2-channel jumps
  void run ();

  // Single iteration with the 2-channel jumps
  VectorXd run_single_iterations (bool verbose = true) const;
};

#endif