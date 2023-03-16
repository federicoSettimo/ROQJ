#include "roqj_gen_qubit.h"

gen_qubit_roqj::gen_qubit_roqj (int N_states, double t_i, double t_f, double dt, int N_copies, bool print_trajectory, int N_traj_print, bool verbose, double threshold) {
  initialize(N_states, t_i, t_f, dt, N_copies, 2, print_trajectory, N_traj_print,verbose,threshold);
  srand(0);
}

Vector2cd gen_qubit_roqj::jump (const Vector2cd &phi1, const Vector2cd &phi2, double lambda1, double lambda2, double z) const {
  if (lambda1 >= -_threshold && lambda2 >= -_threshold) {// Normal jump
    // With probability pjump1, it jumps to the first eigenstate of R
    if (z <= lambda1*_dt)
      return phi1;
    else return phi2;
  }
  else {// Reverse jump ----- Not implemented??
    cerr << "Negative rate - reverse jump. NOT IMPLEMENTED\n";
    cout << lambda1 << ", " << lambda2 << endl;
    exit(EXIT_FAILURE);
  }
}

void gen_qubit_roqj::run () {
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

VectorXd gen_qubit_roqj::run_single_iterations (bool verbose) const {
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

      // The old R without enforcing positivity, we get the eigenstates since they don't change
      MatrixXcd R = J(projector(psi[i]),t) + 0.5*(C(projector(psi[i]), t)*projector(psi[i]) + projector(psi[i])*C(projector(psi[i]), t).transpose());
      ComplexEigenSolver<MatrixXcd> eigs;
      eigs.compute(R);
      VectorXcd eigval = eigs.eigenvalues(), phi, tau, phi_plus;
      MatrixXcd eigvec = eigs.eigenvectors(), Cnew;

      // Computing the new C which will get (hopefully) a positive R
      double lambda1 = real(eigval(0)), lambda2 = real(eigval(1)), lambda, cg = sqrt(real(rho(1,1))), theta = arg(rho(0,1)), beta, norm_inn, lambda_plus;
      complex<double> inner_prod, gamma;
      if (lambda1 < 0.) {
        phi = eigvec.col(0);
        lambda = lambda1;
        phi_plus = eigvec.col(1);
        lambda_plus = lambda2;
      }
      else if (lambda2 < 0.) {
        phi = eigvec.col(1);
        lambda = lambda2;
        phi_plus = eigvec.col(0);
        lambda_plus = lambda2;
      }

      if (lambda1 > 0. && lambda2 > 0.) {
        inner_prod = psi[i].dot(phi);
        norm_inn = norm(inner_prod);
        if (norm_inn == 0) {
          cout << "Warning: orthogonl states\n";
          break;
        }
        beta = lambda/norm_inn;
        gamma = -2.*inner_prod*beta;
        tau = beta*psi[i] + gamma*phi;
        Cnew = tau*((psi[i]).adjoint());
      }
      else Cnew = MatrixXcd::Zero(2,2);

      // Draws a random number and calculates whether the evolution is deterministic or via a jump
      double z = (double)rand()/((double)RAND_MAX);

      if (z < real(R.trace())*_dt){ // Jump
        lambda -= beta*norm_inn;
        lambda_plus += 2.*beta*norm(psi[i].dot(phi_plus))*norm_inn;
        psi[i] = this->jump(phi,phi_plus, lambda,lambda_plus, z);
      }
      else { // Free evolution - now considering also Cnew
        MatrixXcd K = H(t) + 0.5*(C(projector(psi[i]), t).imag() - complex<double>(0.,1.)*(Gamma(t) + C(projector(psi[i]), t).real() ) + Cnew);
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
}