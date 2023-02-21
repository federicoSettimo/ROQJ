#include "roqj_populations.h"

qubit_roqj_pop::qubit_roqj_pop (int N_ensemble, double t_i, double t_f, double dt, int N_copies, bool verbose) {
  srand(time(NULL));
  initialize(N_ensemble, t_i, t_f, dt, N_copies, 2, false, 0, verbose);
  ComplexEigenSolver<MatrixXcd> eigs;
  MatrixXcd R = J(projector(_initial_state),0.5*(_t_f-_t_i)) + 0.5*(C(projector(_initial_state), 0.5*(_t_f-_t_i))*projector(_initial_state) + projector(_initial_state)*C(projector(_initial_state), 0.5*(_t_f-_t_i)).transpose());
  eigs.compute(R);
  MatrixXcd eigvec = eigs.eigenvectors();
  _eig_1 = eigvec.col(0);
  _eig_2 = eigvec.col(1);
  //_eig_1 << 1,0;
  //_eig_2 << 0,1;
  // check that both are not nan
  if (__isnan(real(_eig_2(0))) || __isnan(real(-_eig_2(0)))) {
    if (_eig_1(0) != zero)
      _eig_2 << -conj(_eig_1(1))/conj(_eig_1(0)), 1;
    else _eig_2 << 1, -conj(_eig_1(0))/conj(_eig_1(1));
    _eig_2 = _eig_2.normalized();
  }
  else if (isnan(real(_eig_1(0))) || isnan(real(-_eig_1(0)))) {
    if (_eig_2(0) != zero)
      _eig_1 << -conj(_eig_2(1))/conj(_eig_2(0)), 1;
    else _eig_1 << 1, -conj(_eig_2(0))/conj(_eig_2(1));
    _eig_1 = _eig_1.normalized();
  }
  cout << _eig_1 << endl << _eig_2 << endl;
}

void qubit_roqj_pop::run () {
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

VectorXd qubit_roqj_pop::run_single_iterations (bool verbose) const {
  VectorXd observables(_num_timesteps);
  int n_observable = 0;

  // Allocating _N_states copies of the initial state
  int N_init = _N_states, N_1 = 0, N_2 = 0;

  // Exact solution
  MatrixXcd rho_ex = projector(_initial_state);
  ofstream out_ex, traj;
  if (verbose) {
    out_ex.open("analytic.txt");
    traj.open("trajectories.txt");
  }
  
  // _eig_1 and _eig_2 must be constant, otherwise the ensemble size is not finite
  Vector2cd initial_state_t = _initial_state;
  
  // Time evolution
  for (double t = _t_i; t <= _t_f; t += _dt) {
    // Prints and evolves the exact solution
    if (verbose) {
      out_ex << observable(rho_ex) << endl;
      rho_ex = rho_ex + (-complex<double>(0,1)*comm(H(t),rho_ex) + J(rho_ex,t) - 0.5*anticomm(Gamma(t),rho_ex))*_dt;
      rho_ex /= rho_ex.trace();
    }

    // Average state
    Matrix2cd rho = (double)N_init/((double)_N_states)*projector(initial_state_t) + (double)N_1/((double)_N_states)*projector(_eig_1) + (double)N_2/((double)_N_states)*projector(_eig_2);
    rho = rho/rho.trace();
    observables[n_observable] = observable(rho);
    n_observable++;

    int N_1_old = N_1, N_2_old = N_2, N_init_old = N_init;

    // Cycle on the states in initial_state_t
    MatrixXcd R = J(projector(initial_state_t),t) + 0.5*(C(projector(initial_state_t), t)*projector(initial_state_t) + projector(initial_state_t)*C(projector(initial_state_t), t).transpose());
    ComplexEigenSolver<MatrixXcd> eigs;
    eigs.compute(R);
    VectorXcd eigval = eigs.eigenvalues();
    double lambda_1 = real(eigval[0])*_dt, lambda_2 = real(eigval[1])*_dt;
    MatrixXcd eigvec = eigs.eigenvectors();
    // Lets see if the first eigenvalue is actually _eig_1 or it is eig_2. If eig_2, swap
    if (eigvec.col(0) == _eig_2) {
      double tmp = lambda_1;
      lambda_1 = lambda_2;
      lambda_2 = tmp;
    }
    for (int i = 0; i < N_init; ++i) {
      double z = (double)rand()/((double)RAND_MAX);
      if (z < lambda_1) {
        N_init--; N_1++;
      }
      else if (z < lambda_1 + lambda_2) {
        N_init--; N_2++;
      }
    }

    // Cycle on the states in _eig_1
    R = J(projector(_eig_1),t) + 0.5*(C(projector(_eig_1), t)*projector(_eig_1) + projector(_eig_1)*C(projector(_eig_1), t).transpose());
    eigs.compute(R);
    eigval = eigs.eigenvalues();
    lambda_1 = real(eigval[0])*_dt, lambda_2 = real(eigval[1])*_dt;
    eigvec = eigs.eigenvectors();
    // Lets see if the first eigenvalue is actually _eig_1 or it is eig_2. If eig_2, swap
    if (eigvec.col(0) == _eig_2) {
      double tmp = lambda_1;
      lambda_1 = lambda_2;
      lambda_2 = tmp;
    }
    for (int i = 0; i < N_1_old; ++i) {
      double z = (double)rand()/((double)RAND_MAX);
      if (lambda_1*lambda_2 >= 0) {
        if (z < lambda_2) {
          N_1--; N_2++;
        }
      }
      else if (lambda_1 < 0) {
        if (z <= -lambda_1*(double)N_init_old/((double)N_1_old)) {
          cout << "Reverse jump from eig_1 to initial_state at time " << t << endl;
          N_1--; N_init++;
        }
      }
      else if (lambda_2 < 0) {
        if (z <= -lambda_2*(double)N_2_old/((double)N_1_old)) {
          cout << "Reverse jump from eig_1 to eig_2 at time " << t << endl;
          N_1--; N_2++;
        }
      }
    }

    // Cycle on the states in _eig_2
    R = J(projector(_eig_2),t) + 0.5*(C(projector(_eig_2), t)*projector(_eig_2) + projector(_eig_2)*C(projector(_eig_2), t).transpose());
    eigs.compute(R);
    eigval = eigs.eigenvalues();
    lambda_1 = real(eigval[0])*_dt, lambda_2 = real(eigval[1])*_dt;
    eigvec = eigs.eigenvectors();
    // Lets see if the first eigenvalue is actually _eig_1 or it is eig_2. If eig_2, swap
    if (eigvec.col(0) == _eig_2) {
      double tmp = lambda_1;
      lambda_1 = lambda_2;
      lambda_2 = tmp;
    }
    for (int i = 0; i < N_2_old; ++i) {
      double z = (double)rand()/((double)RAND_MAX);
      if (lambda_1*lambda_2 >= 0) {
        if (z <= lambda_1) {
          N_2--; N_1++;
        }
      }
      else if (lambda_2 < 0) {
        if (z <= -lambda_2*(double)N_init_old/((double)N_2_old)) {
          cout << "Reverse jump from eig_2 to initial_state at time " << t << endl;
          N_2--; N_init++;
        }
      }
      else if (lambda_1 < 0) {
        if (z <= -lambda_1*(double)N_1_old/((double)N_2_old)) {
          cout << "Reverse jump from eig_2 to eig_1 at time " << t << endl;
          N_2--; N_1++;
        }
      }
    }
    MatrixXcd K = H(t) + 0.5*(C(projector(initial_state_t), t).imag() - complex<double>(0.,1.)*(Gamma(t) + C(projector(initial_state_t), t).real() ) );
    initial_state_t -= I*_dt*K*initial_state_t;
    initial_state_t = initial_state_t.normalized();
  }
  return observables;
}