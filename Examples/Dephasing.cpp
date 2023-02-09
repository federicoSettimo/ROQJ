#include "../roqj.h"

using namespace std;
using namespace arma;

const double gamma_z = 5.;

complex<double> I(0,1), one(1,0);
cx_mat sigma_z = {{one,0},{0,-one}}, Id = {{one,0},{0,one}}, sigma_x = {{0,one},{one,0}}, sigma_y = {{0,-I},{I,0}};

double observable (const cx_mat &rho) {
  return real(trace(rho*sigma_z));
}

// Free-evolution effective Hamiltonian
cx_mat H (double t) {
  return cx_mat(2,2,fill::zeros);
}

// J_t(rho)
cx_mat J (const cx_mat &rho, double t) {
  return 0.5*gamma_z*sigma_z*rho*sigma_z;
}

// Gamma(t)
cx_mat Gamma (double t) {
  return .5*gamma_z*Id;
}

// C(t)
cx_mat C (const cx_mat &rho, double t) {
  return (1.-exp(-t))*gamma_z*Id;
}

int main() {
  double tmin = 0., tmax = 10, dt = 0.01;
  int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 15;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, dimH, printTraj, Ntraj);

  jump.set_N_traj_print(Ntraj);

  cx_vec initialState = {sin(M_PI/3.), cos(M_PI/3.)};
  //cx_vec initialState = {1,0};
  //cx_vec initialState = {1./sqrt(2.), 1./sqrt(2.)};
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}