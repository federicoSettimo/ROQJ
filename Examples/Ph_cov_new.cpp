#include "../roqj.h"

using namespace std;
using namespace arma;

complex<double> I(0,1), one(1,0);
cx_mat sigma_p = {{0,one},{0,0}}, sigma_m = {{0,0},{one,0}}, sigma_z = {{one,0}, {0,-one}}, id = {{one,0}, {0,one}}, sigma_x = {{0,one},{one,0}}, sigma_y = {{0,-I},{I,0}};

double gamma_p (double t) {return 1.;}
double gamma_m (double t) {return 1.;}
double gamma_z (double t) {return 0.5*tanh(-t);}
//double gamma_p (double t) {return exp(-2.*t);}
//double gamma_m (double t) {return exp(-t);}
//double gamma_z (double t) {return -0.5*sqrt(gamma_p(t)*gamma_m(t));}
double b (double t) {return 0.5*(1. + erf(2.*sqrt(2.)*(t-1.)));}
//double b (double t) {return 0.;}

cx_mat H (double t) {
  return -.5*b(t)*sigma_z;
}

cx_mat Gamma (double t) {
  return gamma_p(t)*sigma_m*sigma_p + gamma_m(t)*sigma_p*sigma_m + gamma_z(t)*id;
}

/*cx_mat C (const cx_mat &rho, double t) {
  cx_mat P, P_e = {{one,0}, {0,0}}, P_g = {{0,0}, {0,one}};
  double Delta; // maxDelta, aka in the lim z -> \pm 1
  double gp = gamma_p(t), gm = gamma_m(t), gz = gamma_z(t);
  if (gm >= gp) {
    P = P_e;
    Delta = pow(gm-gp,2) + (gm+gp)*(gm+gp-4.*gz) + 2.*(gm-gp)*(gm+gp-2.*gz) + 4.*gz*gz;
  }
  else {
    P = P_g;
    Delta = pow(gm-gp,2) + (gm+gp)*(gm+gp-4.*gz) - 2.*(gm-gp)*(gm+gp-2.*gz) + 4.*gz*gz;
  }
  return sqrt(Delta)*P;
}*/

cx_mat C (const cx_mat &rho, double t) {return .5*(gamma_p(t)+gamma_m(t)+2.*gamma_z(t))*id - 2.*gamma_z(t)*id;}

cx_mat J (const cx_mat &rho, double t) {
  return gamma_p(t)*sigma_p*rho*sigma_m + gamma_m(t)*sigma_m*rho*sigma_p + gamma_z(t)*sigma_z*rho*sigma_z;
}

double observable (const cx_mat &rho) {return real(rho(0,1));}

int main () {
  double tmin = 0., tmax = 10, dt = 0.01;
  int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 5;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, 2, printTraj, Ntraj);

  cx_vec initialState = {sin(M_PI/8.), cos(M_PI/8.)};
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}