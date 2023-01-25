#include "../lib.h"

using namespace std;
using namespace arma;

// Parameters
double gamma1 (double t) {return 1.;}
double gamma2 (double t) {return 1.;}
double gamma3 (double t) {return -.5*tanh(t);}
double gamma_sum (double t) {return gamma1(t)+gamma2(t)+gamma3(t);}
double b (double t) {return 0.5*(1. + erf(2.*sqrt(2.)*(t-1.)));}

double observable (const cx_mat &rho) {
  return real(rho(0,1));
}

complex<double> I(0,1), one(1,0);
cx_mat sigma_z = {{one,0},{0,-one}}, Id = {{one,0},{0,one}}, sigma_x = {{0,one},{one,0}}, sigma_y = {{0,-I},{I,0}};

// Free-evolution effective Hamiltonian
cx_mat H (double t) {
  return -.5*b(t)*sigma_z - Id*I*Gamma(t)*0.5;
}

// J_t(rho)
cx_mat J (const cx_mat &rho, double t) {
  return 0.5*(gamma1(t)*sigma_x*rho*sigma_x + gamma2(t)*sigma_y*rho*sigma_y + gamma3(t)*sigma_z*rho*sigma_z);
}

// Gamma(t)
cx_mat Gamma (double t) {
  return .5*gamma_sum(t)*Id;
}

// C(t)
cx_mat C (const cx_mat &rho, double t) {
  return 0.5*gamma_sum(t)*Id;
}

int main() {
  double tmin = 0., tmax = 10, dt = 0.01;
  int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 5;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj);

  cx_vec initialState = {sin(M_PI/8.), cos(M_PI/8.)};
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}