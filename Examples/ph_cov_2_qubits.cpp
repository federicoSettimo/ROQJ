/*
  Phase covariant dynamics on 2 qubits (non-interacting).
  Works as long as psi1=g,e, but psi2 can differ - why? Where is the simmetry lost?

  Also, for such psi1, it has negative rates also for a CP divisible dynamics (but not when C = 0).
  The problem shoud be in C, only working for the single dephasing

  Even when nothing happens on 1 (gamma_1 = 0) megative rates for psi1 != g,e - how is that possible?
  Also when gamma_z_2 > 0
*/
#include "../roqj.h"

using namespace std;
using namespace Eigen;

double gamma_p_2 (double t) {return 1.;}
double gamma_m_2 (double t) {return 1.;}
double gamma_z_2 (double t) {return -.5*tanh(t);}
//double gamma_p_2 (double t) {return 1.;}
//double gamma_m_2 (double t) {return 1.;}
//double gamma_z_2 (double t) {return -.5*tanh(t);}
double gamma_p_1 (double t) {return 0.;}
double gamma_m_1 (double t) {return 0.;}
double gamma_z_1 (double t) {return 0.;}

//double b (double t) {return 0.5*(1. + erf(2.*sqrt(2.)*(t-1.)));}
double b (double t) {return 0.;}

MatrixXcd H (double t) {
  return -0.5*b(t)*tens(sigma_z, id) - .5*b(t)*tens(id,sigma_z);
}

MatrixXcd J1 (const MatrixXcd &rho, double t) {
  return gamma_p_1(t)*tens(sigma_p, id)*rho*tens(sigma_m, id) + gamma_m_1(t)*tens(sigma_m, id)*rho*tens(sigma_p, id) + gamma_z_1(t)*tens(sigma_z, id)*rho*tens(sigma_z, id);
}

MatrixXcd J2 (const MatrixXcd &rho, double t) {
  return gamma_p_2(t)*tens(id, sigma_p)*rho*tens(id, sigma_m) + gamma_m_2(t)*tens(id, sigma_m)*rho*tens(id, sigma_p) + gamma_z_2(t)*tens(id, sigma_z)*rho*tens(id, sigma_z);
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return J1(rho,t) + J2(rho,t);
  //return gamma_p_1(t)*tens(sigma_p, sigma_p)*rho*tens(sigma_m, sigma_m) + gamma_m_1(t)*tens(sigma_m, sigma_m)*rho*tens(sigma_p, sigma_p) + gamma_z_1(t)*tens(sigma_z, sigma_z)*rho*tens(sigma_z, sigma_z);

}


MatrixXcd Gamma_1 (double t) {
  return gamma_p_1(t)*tens(sigma_m*sigma_p, id) + gamma_m_1(t)*tens(sigma_p*sigma_m, id) + gamma_z_1(t)*tens(id,id);
}

MatrixXcd Gamma_2 (double t) {
  return gamma_p_2(t)*tens(id, sigma_m*sigma_p) + gamma_m_2(t)*tens(id, sigma_p*sigma_m) + gamma_z_2(t)*tens(id,id);
}

MatrixXcd Gamma (double t) {
  return Gamma_1(t) + Gamma_2(t);
  //return gamma_p_1(t)*tens(sigma_m*sigma_p, sigma_m*sigma_p) + gamma_m_1(t)*tens(sigma_p*sigma_m, sigma_p*sigma_m) + gamma_z_1(t)*tens(id,id);
}


MatrixXcd C_1 (const MatrixXcd &rho, double t) {
  double mu, c3, a2 = real(tr_2(rho)(1,1));
  mu = a2 == 0 ? sqrt(gamma_p_1(t)*gamma_m_1(t)) : gamma_m_1(t)*(1.-a2)/a2 - sqrt(gamma_p_1(t)*gamma_m_1(t));
  c3 = gamma_z_1(t) - mu;
  return tens(2.*mu*sigma_p*sigma_m + c3*id, id);
  //return 2.*mu*sigma_p*sigma_m + c3*id;
}

MatrixXcd C_2 (const MatrixXcd &rho, double t) {
  double mu, c3, a2 = real(tr_1(rho)(1,1));
  mu = a2 == 0 ? sqrt(gamma_p_2(t)*gamma_m_2(t)) : gamma_m_2(t)*(1.-a2)/a2 - sqrt(gamma_p_2(t)*gamma_m_2(t));
  c3 = gamma_z_2(t) - mu;
  return tens(id, 2.*mu*sigma_p*sigma_m + c3*id);
  //return 2.*mu*sigma_p*sigma_m + c3*id;
}

MatrixXcd C (const MatrixXcd &rho, double t) {
  return C_1(rho, t) + C_2(rho,t);
  //return tens(C_1(rho,t), C_2(rho,t));
  //return MatrixXcd::Zero(4,4);
}


double observable (const MatrixXcd &rho) {
  return real((sigma_z*tr_1(rho)).trace());
}

int main () {
  double tmin = 0., tmax = 5, dt = 0.01, threshold = 1;
  int N_ensemble = 1000, Ncopies = 2, dimH = 4, Ntraj = 10;
  bool printTraj = true;

  roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, dimH, printTraj, Ntraj, true, threshold);

  //VectorXcd initialState;
  //initialState = VectorXcd::Zero(4);
  //initialState << sin(M_PI/8.), 0, cos(M_PI/8.), 0;
  //initialState << 0., sin(M_PI/8.), 0, cos(M_PI/8.);
  //jump.set_initial_state(initialState);

  Vector2cd psi1;
  psi1 << sin(M_PI/8.), cos(M_PI/8.);
  jump.set_initial_state(tens_state(psi1, excited_state));

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");
  
  return 0;
}