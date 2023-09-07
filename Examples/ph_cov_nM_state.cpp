// Phase covariant violating also P divisibility - no reverse jumps
// Dynamics from Tetteinen - Quantum speed limit and divisibility of the dynamical map
// Tracking with only 1 bit: jumps |psi(t)> -> |0> and |0> -> |psi(t+dt)>
#include "../roqj_state.h"

using namespace std;
using namespace Eigen;

double tmin = 0., tmax = 3., dt = 0.0005, threshold = 1e-3;
vector<Vector2cd> psi_t((int)((tmax-tmin)/dt+1));

double gamma_p (double t) {return exp(-.5*t);}
double gamma_m (double t) {return exp(-.25*t);}
double gamma_z (double t) {return 1.3*exp(-3.*t/8.)*cos(2.*t)*.5;}

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

VectorXcd Phi (const VectorXcd &psi, double t) {
  double a = real(psi(1)), a2 = a*a, gp = gamma_p(t), gm = gamma_m(t), gz = gamma_z(t), sq = sqrt(1.-a2);

  // From |psi>
  if (a2 < 1.) 
    return -(a2*(gp-gz) + gz)/sq*excited_state + (a2*a*(gp-3.*gz) + 3.*a*gz)/(sq*sq)*ground_state;

  // From |0>
  Vector2cd target = psi_t[(int)(t/dt)+1]; // Target = state at the next timestep
  double b = real(target(1)), b2 = b*b, sqb = sqrt(1.-b*b);
  //double b = 0., b2 = b*b, sqb = sqrt(1.-b*b);
  return 2.*b*gp/sqb*excited_state - (gz - b2*(gp+gz))/(sqb*sqb)*ground_state;

  // Just fixed basis
  /*if (a2 > .5)
    return sq*(gm - a2*gm + 3.*a2*gz)/a2*plus_state + (a*(gm-gz) - gm/a)*minus_state;
  return (a2*(gz-gp) - gz)/sq*plus_state + (a2*a*(gp-3.*gz) + 3.*a*gz)/(sq*sq)*minus_state;*/
  /*if (a2 > .5)
    return ((sqrt(1 - pow(a,2))*(gm - pow(a,2)*gm + 3*pow(a,2)*gz))/pow(a,2))*excited_state
      + (-(gm/a) + a*(gm - gz))*ground_state;
  return ((-gz + pow(a,2)*(-gp + gz))/sqrt(1 - pow(a,2)))*excited_state 
    + (-((pow(a,3)*(gp - 3*gz) + 3*a*gz)/(-1 + pow(a,2))))*ground_state;*/
}

//double observable (const MatrixXcd &rho) {return abs(rho(0,1));}
double observable (const MatrixXcd &rho) {return real((rho*sigma_x).trace());}

int main () {
  int N_ensemble = 1000, Ncopies = 3, dimH = 2, Ntraj = 10;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj, true, threshold);
  Vector2cd initialState;
  initialState << 0.8, 0.7;
  jump.set_initial_state(initialState);

  // First thing: generating the deterministic evolution
  int i = 1;
  psi_t[0] = initialState.normalized();
  for (double t = tmin; t < tmax; t += dt) {
    Matrix2cd K = H(t) - .5*I*Gamma(t);
    psi_t[i] = psi_t[i-1] - I*dt*K*psi_t[i-1] - .5*dt*Phi(psi_t[i-1],t);
    psi_t[i].normalize();
    i++;
  }

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}