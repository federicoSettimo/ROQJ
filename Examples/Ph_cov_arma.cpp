#include "../roqj_arma.h"
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace arma;

complex<double> I(0,1), one(1,0);
cx_mat sigma_p = {{0,one},{0,0}}, sigma_m = {{0,0},{one,0}}, sigma_z = {{one,0}, {0,-one}}, id = {{one,0}, {0,one}}, sigma_x = {{0,one},{one,0}}, sigma_y = {{0,-I},{I,0}};
double gamma_p (double t) {return 1.;}
double gamma_m (double t) {return 1.;}
double gamma_z (double t) {return -.5*tanh(t);}
//double gamma_p (double t) {return exp(-2.*t);}
//double gamma_m (double t) {return exp(-t);}
//double gamma_z (double t) {return -0.5*sqrt(gamma_p(t)*gamma_m(t));}
double b (double t) {return 0.5*(1. + erf(2.*sqrt(2.)*(t-1.)));}
//double b (double t) {return 0.;}

cx_mat H (double t) {
  return -0.5*b(t)*sigma_z;
}

cx_mat J (const cx_mat &rho, double t) {
  return gamma_p(t)*sigma_p*rho*sigma_m + gamma_m(t)*sigma_m*rho*sigma_p + gamma_z(t)*sigma_z*rho*sigma_z;
}

cx_mat Gamma (double t) {
  return gamma_p(t)*sigma_m*sigma_p + gamma_m(t)*sigma_p*sigma_m + gamma_z(t)*id;
}

cx_mat C (const cx_mat &rho, double t) {
  double mu, c3, a2 = real(rho.at(1,1));
  mu = a2 == 0 ? sqrt(gamma_p(t)*gamma_m(t)) : gamma_m(t)*(1.-a2)/a2 - sqrt(gamma_p(t)*gamma_m(t));
  c3 = gamma_z(t) - mu;
  return 2.*mu*sigma_p*sigma_m + c3*id;
}

// C for R1 for the Eternally non-Markov
//cx_mat C (const cx_mat &rho, double t) {return .5*(gamma_p(t)+gamma_m(t)+2.*gamma_z(t))*id;}

//double observable (const cx_mat &rho) {return real(rho(0,1));}
double observable (const cx_mat &rho) {return real(trace(rho*sigma_x));}

int main () {
  double tmin = 0., tmax = 5, dt = 0.01;
  int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 10;
  bool printTraj = true;

  roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, 2, printTraj, Ntraj);
  jump.set_N_traj_print(Ntraj);

  cx_vec initialState = {sin(M_PI/8.), cos(M_PI/8.)};
  jump.set_initial_state(initialState);

  auto start = high_resolution_clock::now();
  jump.run();
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<seconds>(stop - start);
  cout << "\tExecution time: " << duration.count() << " s\n";

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}