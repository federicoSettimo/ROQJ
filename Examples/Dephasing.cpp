#include "../roqj.h"

using namespace std;
using namespace Eigen;

const double gamma_z = 5.;
double observable (const MatrixXcd &rho) {
  return real((rho*sigma_z).trace());
}

// Free-evolution effective Hamiltonian
MatrixXcd H (double t) {
  return MatrixXcd::Zero(2,2);
}

// J_t(rho)
MatrixXcd J (const MatrixXcd &rho, double t) {
  return 0.5*gamma_z*sigma_z*rho*sigma_z;
}

// Gamma(t)
MatrixXcd Gamma (double t) {
  return .5*gamma_z*id;
}

// C(t)
MatrixXcd C (const MatrixXcd &rho, double t) {
  return (1.-exp(-t))*gamma_z*id;
}

int main() {
  double tmin = 0., tmax = 10, dt = 0.01;
  int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 15;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, dimH, printTraj, Ntraj);

  jump.set_N_traj_print(Ntraj);

  VectorXcd initialState(2);
  initialState << sin(M_PI/3.), cos(M_PI/3.);
  //initialState << 1,0;
  //initialState << 1./sqrt(2.), 1./sqrt(2.);
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}