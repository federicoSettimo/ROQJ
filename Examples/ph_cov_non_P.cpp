// Non-P divisible phase covariant. Model from Tetteinen - Quantum Speed Limit and Divisibility of the Dynamical Map
#include "../roqj_pop.h"

using namespace std;
using namespace Eigen;

double gamma_p (double t) {return exp(-.5*t);}
double gamma_m (double t) {return exp(-.25*t);}
double gamma_z (double t) {return 2.*0.5*exp(-3.*t/8.)*cos(2.*t);}
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

MatrixXcd C (const VectorXcd &psi, double t) {
  double mu, c3, a2 = norm(psi(1));
  mu = a2 == 0. ? sqrt(gamma_p(t)*gamma_m(t)) : gamma_m(t)*(1.-a2)/a2 - sqrt(gamma_p(t)*gamma_m(t));
  //mu = a2 == 1. ? -sqrt(gamma_p(t)*gamma_m(t)) : -gamma_p(t)*a2/(1.-a2) + sqrt(gamma_p(t)*gamma_m(t));
  //if (a2 == 0.) mu = sqrt(gamma_p(t)*gamma_m(t));
  //else if (a2 == 1.) mu = -sqrt(gamma_p(t)*gamma_m(t));
  //else mu = .5*(gamma_m(t)*(1.-a2)/a2 - sqrt(gamma_p(t)*gamma_m(t))) + .5*(-gamma_p(t)*a2/(1.-a2) + sqrt(gamma_p(t)*gamma_m(t)));
  c3 = gamma_z(t) - mu;
  return 2.*mu*sigma_p*sigma_m + c3*id;
}

double observable (const MatrixXcd &rho) {return real((rho*sigma_z).trace());}

int main () {
  double tmin = 0., tmax = 5, dt = 0.01;
  int N_ensemble = 10000, Ncopies = 5, dimH = 2, Ntraj = 10;
  bool printTraj = true;

  qubit_roqj_pop jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj);

  Vector2cd initialState;
  initialState << sin(M_PI/8.), cos(M_PI/8.);
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");
  jump.get_trajectories();

  return 0;
}