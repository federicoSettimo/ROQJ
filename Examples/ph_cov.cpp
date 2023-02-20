#include "../roqj.h"
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace Eigen;

double gamma_p (double t) {return 1.;}
double gamma_m (double t) {return 1.;}
double gamma_z (double t) {return -.5*tanh(t);}
//double gamma_p (double t) {return exp(-2.*t);}
//double gamma_m (double t) {return exp(-t);}
//double gamma_z (double t) {return -0.5*sqrt(gamma_p(t)*gamma_m(t));}
//double gamma_z (double t) {return -sqrt(gamma_p(t)*gamma_m(t));}
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
  double mu, c3, a2 = real(rho(1,1));
  mu = a2 == 0 ? sqrt(gamma_p(t)*gamma_m(t)) : gamma_m(t)*(1.-a2)/a2 - sqrt(gamma_p(t)*gamma_m(t));
  c3 = gamma_z(t) - mu;
  return 2.*mu*sigma_p*sigma_m + c3*id;
}

// C for R1 for the Eternally non-Markov
//MatrixXcd C (const MatrixXcd &rho, double t) {return .5*(gamma_p(t)+gamma_m(t)+2.*gamma_z(t))*id;}

double observable (const MatrixXcd &rho) {return real(rho(0,1));}
//double observable (const MatrixXcd &rho) {return real((rho*sigma_z).trace());}

int main () {
  double tmin = 0., tmax = 5, dt = 0.01;
  int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 10;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj);

  Vector2cd initialState;
  initialState << sin(M_PI/8.), cos(M_PI/8.);
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}