/*
  Phase covariant dynamics on 2 qubits (non-interacting).

  Fixed J = (J1 tens id) + (id tens J2)

  Works only for factorised initial states and only one qubit evolving. If both evolve, jumps are to entangled states and
  therefore the positivity of R is no longer guaranteed.
*/
#include "../roqj.h"
#include <Unsupported/Eigen/KroneckerProduct>

using namespace std;
using namespace Eigen;

double gamma_p_1 (double t) {return 1.;}
double gamma_m_1 (double t) {return 1.;}
double gamma_z_1 (double t) {return -.5*tanh(t);}
double gamma_p_2 (double t) {return gamma_p_1(t);}
double gamma_m_2 (double t) {return gamma_m_1(t);}
double gamma_z_2 (double t) {return gamma_z_1(t);}
//double gamma_p_2 (double t) {return 0.;}
//double gamma_m_2 (double t) {return 0.;}
//double gamma_z_2 (double t) {return 0.;}

//double b (double t) {return 0.5*(1. + erf(2.*sqrt(2.)*(t-1.)));}
double b (double t) {return 0.;}

MatrixXcd H (double t) {
  //return -0.5*b(t)*tens(sigma_z, id) - .5*b(t)*tens(id,sigma_z);
  return MatrixXcd::Zero(4,4);
}

MatrixXcd J1 (const MatrixXcd &rho, double t) {
  return gamma_p_1(t)*tens(sigma_p, id)*rho*tens(sigma_m, id) + gamma_m_1(t)*tens(sigma_m, id)*rho*tens(sigma_p, id) + gamma_z_1(t)*tens(sigma_z, id)*rho*tens(sigma_z, id);
}

MatrixXcd J2 (const MatrixXcd &rho, double t) {
  return gamma_p_2(t)*tens(id, sigma_p)*rho*tens(id, sigma_m) + gamma_m_2(t)*tens(id, sigma_m)*rho*tens(id, sigma_p) + gamma_z_2(t)*tens(id, sigma_z)*rho*tens(id, sigma_z);
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return J1(rho,t) + J2(rho,t);
}


MatrixXcd Gamma_1 (double t) {
  return gamma_p_1(t)*tens(sigma_m*sigma_p, id) + gamma_m_1(t)*tens(sigma_p*sigma_m, id) + gamma_z_1(t)*tens(id,id);
}

MatrixXcd Gamma_2 (double t) {
  return gamma_p_2(t)*tens(id, sigma_m*sigma_p) + gamma_m_2(t)*tens(id, sigma_p*sigma_m) + gamma_z_2(t)*tens(id,id);
}

MatrixXcd Gamma (double t) {
  return Gamma_1(t) + Gamma_2(t);
}


MatrixXcd C_1 (const MatrixXcd &rho, double t) {
  //return .5*(gamma_m_1(t)+gamma_p_1(t)+ gamma_z_1(t))*id; // C from the paper

  double mu, c3, a2 = real(tr_2(rho)(1,1));
  mu = abs(a2) <= 1e-4 ? sqrt(gamma_p_1(t)*gamma_m_1(t)) : gamma_m_1(t)*(1.-a2)/a2 - sqrt(gamma_p_1(t)*gamma_m_1(t));
  c3 = gamma_z_1(t) - mu;
  return 2.*mu*sigma_p*sigma_m + c3*id;
}

MatrixXcd C_2 (const MatrixXcd &rho, double t) {
  //return .5*(gamma_m_2(t)+gamma_p_2(t)+gamma_z_2(t))*id; // C from the paper

  double mu, c3, a2 = real(tr_1(rho)(1,1));
  mu = abs(a2) <= 1e-4 ? sqrt(gamma_p_2(t)*gamma_m_2(t)) : gamma_m_2(t)*(1.-a2)/a2 - sqrt(gamma_p_2(t)*gamma_m_2(t));
  c3 = gamma_z_2(t) - mu;
  return 2.*mu*sigma_p*sigma_m + c3*id;
}

MatrixXcd C (const MatrixXcd &rho, double t) {
  return tens(C_1(rho,t), id) + tens(id, C_2(rho,t));
}


double observable (const MatrixXcd &rho) {
  return real((sigma_z*tr_2(rho)).trace());
}


int main () {
  double tmin = 0., tmax = 5., dt = 0.01, threshold = 1e-5;
  int N_ensemble = 1000, Ncopies = 3, dimH = 4, Ntraj = 10;
  bool printTraj = true;

  roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, dimH, printTraj, Ntraj, true, threshold);

  VectorXcd entangled(4);
  entangled << 1./sqrt(2.), 0., 0., 1./sqrt(2.);
  jump.set_initial_state(entangled);

  Vector2cd psi;
  psi << sin(M_PI/8.), cos(M_PI/8.);
  //VectorXcd factorized = tens_state(plus_state,plus_state);
  VectorXcd factorized = tens_state(psi,psi);

  //jump.set_initial_state(entangled);
  jump.set_initial_state(factorized);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");
  
  return 0;
}