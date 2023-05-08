/*
  Phase covariant dynamics on 2 qubits, each with separate environment (J = J1 tens id + id tens J2).
  Interaction generating entanglement between the two, observing the entanglement of formation for the trajectories
*/
#include "../roqj.h"
#include <Unsupported/Eigen/KroneckerProduct>

using namespace std;
using namespace Eigen;

///*
double gamma_p (double t) {return exp(-.25*t);}
double gamma_m (double t) {return exp(-.25*t);}
double gamma_z (double t) {return 0.;}
//*/
/*
double gamma_p (double t) {return 1.;}
double gamma_m (double t) {return 1.;}
double gamma_z (double t) {return -.5*tanh(t);}
//*/

MatrixXcd H (double t) {
  return tens(sigma_p,sigma_m) + tens(sigma_m,sigma_p);
  //return MatrixXcd::Zero(4,4);
}


MatrixXcd J (const MatrixXcd &rho, double t) {
  return gamma_p(t)*(tens(sigma_p,id)*rho*tens(sigma_m,id) + tens(id,sigma_p)*rho*tens(id,sigma_m)) +
         gamma_m(t)*(tens(sigma_m,id)*rho*tens(sigma_p,id) + tens(id,sigma_m)*rho*tens(id,sigma_p)) +
         gamma_z(t)*(tens(sigma_z,id)*rho*tens(sigma_z,id) + tens(id,sigma_z)*rho*tens(id,sigma_z));
}


MatrixXcd Gamma (double t) {
  return gamma_p(t)*tens(sigma_m*sigma_p, id) + gamma_m(t)*tens(sigma_p*sigma_m, id) + gamma_z(t)*tens(id,id)
       + gamma_p(t)*tens(id, sigma_m*sigma_p) + gamma_m(t)*tens(id, sigma_p*sigma_m) + gamma_z(t)*tens(id,id);
}

MatrixXcd C (const VectorXcd &psi, double t) {
  Matrix4cd Ppsi = projector(psi), JJ = J(Ppsi,t);
  complex<double> J11 = psi.dot(JJ*psi);
  return -J11*Ppsi + 2.*J11*Ppsi - JJ*Ppsi; // The correct C: ok for P div, not ok for enm (mathematica confirms a <0 eigenvalue)
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