/*
  Phase covariant dynamics on 2 qubits, each with separate environment (J = J1 tens id + id tens J2).
  Interaction generating entanglement between the two, observing the entanglement of formation for the trajectories
*/
#include "../roqj.h"
#include <Unsupported/Eigen/KroneckerProduct>

using namespace std;
using namespace Eigen;

/*
double gamma_p (double t) {return exp(-.25*t);}
double gamma_m (double t) {return exp(-.25*t);}
double gamma_z (double t) {return 0.;}
//*/
///*
double gamma_p (double t) {return 1.;}
double gamma_m (double t) {return 1.;}
double gamma_z (double t) {return -.5*tanh(t);}
//*/

MatrixXcd H (double t) {
  return tens(sigma_p,sigma_m) + tens(sigma_m,sigma_p);
  //return MatrixXcd::Zero(4,4);
}


Matrix2cd J1 (const Matrix2cd &rho, double t) {
  return gamma_p(t)*sigma_p*rho*sigma_m + gamma_m(t)*sigma_m*rho*sigma_p + gamma_z(t)*sigma_z*rho*sigma_z;
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return tens(J1(tr_2(rho),t),id) + tens(id,J1(tr_1(rho),t));
}


MatrixXcd Gamma (double t) {
  return gamma_p(t)*tens(sigma_m*sigma_p, id) + gamma_m(t)*tens(sigma_p*sigma_m, id) + gamma_z(t)*tens(id,id)
       + gamma_p(t)*tens(id, sigma_m*sigma_p) + gamma_m(t)*tens(id, sigma_p*sigma_m) + gamma_z(t)*tens(id,id);
}

MatrixXcd C (const VectorXcd &psi, double t) {
  Matrix4cd Ppsi = projector(psi), JJ = J(Ppsi,t);
  complex<double> J11 = psi.dot(JJ*psi);
  //return -J11*Ppsi + 2.*J11*Ppsi - JJ*Ppsi; // The correct C: ok for P div, not ok for enm (mathematica confirms a <0 eigenvalue)

  return -real(J11)*Ppsi + 2.*Ppsi*JJ*Ppsi - 2.*JJ*Ppsi - 2.*Ppsi*JJ; // With this extra term, ok for enm, not ok for P div dynamics
  // Also, now the jumps preserve the entropy
}


double observable (const MatrixXcd &rho) {
  return entropy(tr_1(rho));
  //return real((sigma_x*tr_1(rho)).trace());
}


int main () {
  double tmin = 0., tmax = 5., dt = 0.01, threshold = 1e-5;
  int N_ensemble = 100, Ncopies = 1, dimH = 4, Ntraj = 10;
  bool printTraj = true;

  roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, dimH, printTraj, Ntraj, true, threshold);

  Vector4cd psi = tens_state(plus_state,plus_state);
  jump.set_initial_state(psi);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");
  
  return 0;
}