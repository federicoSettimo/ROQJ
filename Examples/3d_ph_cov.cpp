/*
  3d phase covariant

  At the moment, it works fine for the coherences but not for the populations
  Also, only some initial states are ok
*/
#include "../roqj.h"

using namespace std;
using namespace Eigen;

double gamma_deph (double t) {return -tanh(t)/3.;}

MatrixXcd H (double t) {
  return MatrixXcd::Zero(3,3);
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  double psi0 = real(rho(0,0)), psi1 = real(rho(1,1)), psi2 = real(rho(2,2)), g = gamma_deph(t);
  MatrixXcd JJ = -g*rho;
  JJ(0,0) = g*psi0 + psi1 + psi2;
  JJ(1,1) = psi0 + g*psi1 + psi2;
  JJ(2,2) = psi0 + psi1 + g*psi2;
  return JJ;
}

MatrixXcd Gamma (double t) {
  return 2.*(1. + gamma_deph(t))*MatrixXcd::Identity(3,3);
}

MatrixXcd C (const VectorXcd &psi, double t) {
  double g = gamma_deph(t), psi02 = norm(psi(0)), psi12 = norm(psi(1)), psi22 = norm(psi(2)), c2;
  if (psi02 != 0.)
    c2 = 3.*g + psi12/psi02 + psi22/psi02; // c2ub
  else c2 = (g + 2.*g*psi22*psi22 + 3.*(psi22 + psi22*psi22))/(-1. - 2.*psi22*psi22 + 4.*g*psi22*(psi22-1.));//c2lb if psi02 = 0

  MatrixXcd CC = MatrixXcd::Zero(3,3);

  CC(2,2) = c2;
  CC(1,1) = c2;
  CC(0,0) = 2.*g - c2;

  return CC;
}

double observable (const MatrixXcd &rho) {
  return real(rho(0,1));
}

int main () {
  double tmin = 0., tmax = 5, dt = 0.01, threshold = 1e-3;
  int N_ensemble = 1000, Ncopies = 3, dimH = 3, Ntraj = 10;
  bool printTraj = true;

  roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, dimH, printTraj, Ntraj, true, threshold);

  Vector3cd initialState;
  initialState << 1., 1., 1.;
  initialState.normalize();
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}