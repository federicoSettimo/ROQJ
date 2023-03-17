/*
  Using a RO that should work for any P divisible qubit dynamics, as long as the pre-jump state |psi>
  and the negative eigenstate of J are not orthogonal.
  It is however very bad: at each timestep and for each state it calculates the eigens of J,
  gets the new R and calculates eigs again - not good

  Note: it can be improved with a new class of roqj:
  Some C just to modifies the eigenstates (e.g. jumps to e,g without enforcing positivity),
  get the new R (also the old one works fine), and does the modifications here done in C
  only in the jump method: eigenstates are calculted only once and not modified by this C
*/

#include "../roqj.h"
#include "../roqj_gen_qubit.h"

using namespace std;
using namespace Eigen;

double gamma_p (double t) {return 1.;}
double gamma_m (double t) {return 1.;}
double gamma_z (double t) {return -.5*tanh(t);}

//double b (double t) {return 0.5*(1. + erf(2.*sqrt(2.)*(t-1.)));}
double b (double t) {return 0.;}

MatrixXcd H (double t) {
  return -0.5*b(t)*sigma_z;
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return gamma_p(t)*sigma_p*rho*sigma_m + gamma_m(t)*sigma_m*rho*sigma_p + gamma_z(t)*sigma_z*rho*sigma_z;
}

MatrixXcd Gamma (double t) {
  return gamma_p(t)*sigma_m*sigma_p + gamma_m(t)*sigma_p*sigma_m + gamma_z(t)*id;
}

MatrixXcd C (const MatrixXcd &rho, double t) {
  // Using the optimized version - working on z
  //return MatrixXcd::Zero(2,2);

  // Using the optimized version, jumps to |g,e> - not working
  double a2 = abs(rho(1,1)), phi = arg(rho(0,1)), r = abs(J(rho,t)(0,1)), theta = arg(J(rho,t)(0,1)), sphi = sin(phi);
  if (a2 == 0 || a2 == 1)
    return MatrixXcd::Zero(2,2);
  MatrixXcd C = MatrixXcd::Zero(2,2);
  if (sphi == 0) {
    C(0,0) = -2.*r/(sqrt(a2)*sqrt(1.-a2));
    C(1,0) = 2.*r*sin(theta)/(1.-a2);
  }
  else {
    C(0,0) = -2.*r*sin(theta)/(sqrt(a2)*sqrt(1.-a2)*sphi);
    C(1,0) = 2.*r*sin(theta-phi)/((1.-a2)*sphi);
  }
  return C;

  // Using the unoptimized version - working
  /*ComplexEigenSolver<MatrixXcd> eigs;
  eigs.compute(J(rho,t));
  VectorXcd eigval = eigs.eigenvalues(), phi, psi, tau;
  MatrixXcd eigvec = eigs.eigenvectors();

  double lambda1 = real(eigval(0)), lambda2 = real(eigval(1)), lambda, cg = sqrt(real(rho(1,1))), theta = arg(rho(0,1)), beta, norm_inn;
  complex<double> inner_prod, gamma;
  if (lambda1 < 0.) {
    phi = eigvec.col(0);
    lambda = lambda1;
  }
  else if (lambda2 < 0.) {
    phi = eigvec.col(1);
    lambda = lambda2;
  }
  else return MatrixXcd::Zero(2,2);

  psi = cg*ground_state + sqrt(1.-cg*cg)*exp(I*theta)*excited_state;
  inner_prod = psi.dot(phi);
  norm_inn = norm(inner_prod);
  if (norm_inn == 0) {
    cout << "Wrong: orthogonl states\n";
    return MatrixXcd::Zero(2,2);
  }
  beta = lambda/norm_inn;
  gamma = -2.*inner_prod*beta;
  tau = beta*psi + gamma*phi;
  return tau*psi.adjoint();*/
}

//double observable (const MatrixXcd &rho) {return real(rho(0,1));}
double observable (const MatrixXcd &rho) {return real((rho*sigma_z).trace());}

int main () {
  double tmin = 0., tmax = 5, dt = 0.01, threshold = 1e-5;
  int N_ensemble = 1000, Ncopies = 3, dimH = 2, Ntraj = 10;
  bool printTraj = true;

  //qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj);
  gen_qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj, true, threshold);

  Vector2cd initialState;
  initialState << sin(M_PI/8.), cos(M_PI/8.);
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}