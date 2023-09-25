// Calculating det. evol for some particular inital state
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>
#include <vector>
#include <cstring>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace Eigen;

static complex<double> I(0,1), one(1,0), zero(0,0);
static Eigen::MatrixXcd sigma_x {{0,1},{1,0}};
static Eigen::MatrixXcd sigma_y {{0,-I},{I,0}};
static Eigen::MatrixXcd sigma_z {{1,0},{0,-1}};
static Eigen::MatrixXcd sigma_p {{0,1},{0,0}};
static Eigen::MatrixXcd sigma_m {{0,0},{1,0}};
static Eigen::MatrixXcd id {{1,0},{0,1}};





// ------------------------- SOME STATES -------------------------
static Eigen::VectorXcd ground_state {{0.,1.}};
static Eigen::VectorXcd excited_state {{1.,0.}};
static Eigen::VectorXcd plus_state {{1./sqrt(2.),1./sqrt(2.)}};
static Eigen::VectorXcd minus_state {{1./sqrt(2.),-1./sqrt(2.)}};

double tmin = 0., tmax = 100., dt = 0.5, threshold = 1e-3;
vector<Vector2cd> psi_t((int)((tmax-tmin)/dt+1));

double beta (double t) {return 1.;}
double gamma_m (double t) {return 1.;}

MatrixXcd H (double t) {
  return beta(t)*sigma_y;
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return gamma_m(t)*sigma_m*rho*sigma_p;
}

MatrixXcd Gamma (double t) {
  return gamma_m(t)*sigma_m;
}

VectorXcd Phi (const VectorXcd &psi, double t) {
  double b = beta(t), gm = gamma_m(t), mu = real(psi(1));

  if (abs(mu) == 1.) 
    return (-2*b*dt)*excited_state + mu*ground_state;

  if (abs(mu) == 0.)
    return 2.*excited_state + (-2*b*dt*(-2 + gm))*ground_state;

  return (sqrt(1 - pow(mu,2)))*excited_state + (mu - 2*gm*mu + 2*gm*pow(mu,3))*ground_state;
}

double observable (const MatrixXcd &rho) {return real((rho*sigma_z).trace());}

Matrix2cd projector (const Vector2cd &psi) {return psi*psi.adjoint();}

int main () {
  ofstream out;
  out.open("det_evol.txt");
  out << tmin << endl << tmax << endl << dt << endl;

  Vector2cd psi;
  psi << 0., 1.;
  psi.normalize();

  for (double t = tmin; t <= tmax; t += dt) {
    Vector2cd phi = Phi(psi,t);
    out << observable(projector(psi)) << endl;
    out << observable(projector(phi.normalized())) << endl;
    out << phi.norm() << endl;
    out << observable(projector(phi)) << endl;
    Matrix2cd K = H(t) - .5*I*Gamma(t);
    psi -= I*dt*K*psi + .5*dt*phi;
    psi.normalize();
  }
  
  return 0;
}