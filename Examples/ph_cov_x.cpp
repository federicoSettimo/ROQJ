/*
  Again the generic phase covariant, but now supposing that the driving is sigma_x instead of sigma_z.
  Thus using a C such that the eigenstates of R are the eigenstates of sigma_x |+,->.
  Writing C = sum_{i=0}^3 a_i sigma_i, and a_i s.t. R = c1 id + c2 sigma_x and eigval(R)>0
*/
#include "../roqj.h"

using namespace std;
using namespace Eigen;

//double gamma_p (double t) {return 1.;}
//double gamma_m (double t) {return 1.;}
//double gamma_z (double t) {return -.5*tanh(t);}
double gamma_p (double t) {return exp(-t);}
double gamma_m (double t) {return exp(-t);}
double gamma_z (double t) {return -0.5*sqrt(gamma_p(t)*gamma_m(t));}
//double gamma_z (double t) {return -sqrt(gamma_p(t)*gamma_m(t));}
double b (double t) {return 2.*0.5*(1. + erf(2.*sqrt(2.)*(t-1.)));}
//double b (double t) {return 0.;}

MatrixXcd H (double t) {
  return -0.5*b(t)*sigma_x;
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return gamma_p(t)*sigma_p*rho*sigma_m + gamma_m(t)*sigma_m*rho*sigma_p + gamma_z(t)*sigma_z*rho*sigma_z;
}

MatrixXcd Gamma (double t) {
  return gamma_p(t)*sigma_m*sigma_p + gamma_m(t)*sigma_p*sigma_m + gamma_z(t)*id;
}

MatrixXcd C (const MatrixXcd &rho, double t) {
  double z, a = sqrt(real(rho(1,1))), a2 = a*a, gm = gamma_m(t), gp = gamma_p(t), gz = gamma_z(t);
  z = gm*(1-a2) - gp*a2 + gz*(2.*a2-1);
  return z*sigma_z;
}

double observable (const MatrixXcd &rho) {return real((rho*sigma_y).trace());}

int main () {
  double tmin = 0., tmax = 5, dt = 0.01;
  int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 10;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj);

  Vector2cd initialState;
  initialState << sin(M_PI/8.), cos(M_PI/8.);
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average_2.txt");
  jump.get_error_observable("error_2.txt");

  return 0;
}