/*
  3d phase covariant

  Something like the eternally non-Markovian, with gamma = - tanh(t)/3.
  Unlike qubit, the coherences tend to zero
*/
#include "../roqj.h"

using namespace std;
using namespace Eigen;

double gamma_deph (double t) {return -tanh(t)/3.;}

MatrixXcd H (double t) {
  return MatrixXcd::Zero(3,3);
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  complex<double> psi0 = rho(0,0), psi1 = rho(1,1), psi2 = rho(2,2);
  double g = gamma_deph(t);
  MatrixXcd JJ = -g*rho;
  JJ(0,0) = 2.*g*psi0 + psi1 + psi2;
  JJ(1,1) = psi0 + 2.*g*psi1 + psi2;
  JJ(2,2) = psi0 + psi1 + 2.*g*psi2;
  return JJ;
}

MatrixXcd Gamma (double t) {
  return (1. + 2.*gamma_deph(t))*MatrixXcd::Identity(3,3);
}

MatrixXcd C (const VectorXcd &psi, double t) {
  double g = gamma_deph(t), psi0 = norm(psi(0)), psi1 = norm(psi(1)), psi2 = norm(psi(2)), c;
  
  if (psi0 != 0.)
    c = 4.*g + (psi1 + psi2)/psi0;
  else if (psi1*psi1 + 6.*g*psi1*psi2 + psi2*psi2 != 0.)
    c = (-2.*g*(psi1*psi1+psi2*psi2) - (1.+3.*g)*psi1*psi2)/(psi1*psi1 + 6.*g*psi1*psi2 + psi2*psi2);
  else
    c = 0.;

  MatrixXcd CC = MatrixXcd::Zero(3,3);

  CC(2,2) = c;
  CC(1,1) = c;
  CC(0,0) = 2.*g - c;

  return CC;
}


double observable (const MatrixXcd &rho) {
  //return real(rho(0,1));
  return real(rho(1,1) - rho(2,2));
}

int main () {
  double tmin = 0., tmax = 5, dt = 0.01, threshold = 1e-3;
  int N_ensemble = 10000, Ncopies = 3, dimH = 3, Ntraj = 10;
  bool printTraj = true;

  roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, dimH, printTraj, Ntraj, true, threshold);

  Vector3cd initialState;
  initialState << 3., 2., 1.;
  initialState.normalize();
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}