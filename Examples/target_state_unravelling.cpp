// Doing the unravelling using |Phi><psi| instead of C
// Plotting det evolution and evolution of Phi
// Using as a model gamma_pm<0, fixed basis |+,-> for jumps
// Since there is a freedom for Phi, plotting the 2 limit cases and the inbetween one

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

// ------------------------- PAULI MATRICES -------------------------
static complex<double> I(0,1), one(1,0), zero(0,0);
static Eigen::MatrixXcd sigma_x {{0,1},{1,0}};
static Eigen::MatrixXcd sigma_y {{0,-I},{I,0}};
static Eigen::MatrixXcd sigma_z {{1,0},{0,-1}};
static Eigen::MatrixXcd sigma_p {{0,1},{0,0}};
static Eigen::MatrixXcd sigma_m {{0,0},{1,0}};
static Eigen::MatrixXcd id {{1,0},{0,1}};
Matrix2cd sigma_p_pm {{1,-1},{1,-1}}, sigma_m_pm {{1,1},{-1,-1}}, sigma_z_pm {{0,1},{1,0}}, sigma_x_pm {{1,0},{0,-1}}, sigma_y_pm {{0,I},{-I,0}};

// ------------------------- SOME STATES -------------------------
static Eigen::VectorXcd plus_state {{0.,1.}};
static Eigen::VectorXcd minus_state {{1.,0.}};

// Note: it must be gp + gm + gz >= 0 at all times
double kappa = .35;
double gamma_p (double t) {return exp(-.1*t)*(kappa + (1.-kappa)*exp(-.25*t)*cos(2.*t));}
double gamma_m (double t) {return gamma_p(t);}
double gamma_z (double t) {return 1.;}

Matrix2cd H (double t) {
  return MatrixXcd::Zero(2,2);
}

Matrix2cd J (const MatrixXcd &rho, double t) {
  return gamma_p(t)*sigma_p_pm*rho*sigma_m_pm + gamma_m(t)*sigma_m_pm*rho*sigma_p_pm + gamma_z(t)*sigma_z_pm*rho*sigma_z_pm;
}

Matrix2cd Gamma (double t) {
  return gamma_p(t)*sigma_m_pm*sigma_p_pm + gamma_m(t)*sigma_p_pm*sigma_m_pm + gamma_z(t)*id;
}

Vector2cd Phi (const VectorXcd &psi, double t) {
  double a = real(psi(1)), cp, cm; // a = <psi|->
  double a2 = a*a, g = gamma_m(t), gz = gamma_z(t), sq = sqrt(1.-a2);

  if (a2 > .9) {
    cp = sq*(2.*g + 8.*a2*g + gz - 3.*a2*gz)/(2.*a2);;
    cm = -(2.*g + gz - a2*gz)/(2.*a);
  }
  else {
    cp = (-2.*g - a2*gz)/(2.*sq);
    cm = (10.*a*g - 8.*a2*a*g - 2.*a*gz + 3.*a2*gz) /(2.*sq*sq);
  }

  return cp*plus_state + cm*minus_state;
}

Matrix2cd projector (const Vector2cd &psi) {return psi*psi.adjoint();}

int main () {
  double tmin = 0., tmax = 2., dt = .01;
  int npoints = (int)((tmax-tmin)/dt+1);

  Vector2cd initialState;
  initialState << 0.2, -0.7;

  ofstream out;
  out.open("target_state_unravelling.txt");
  out << tmax << endl << dt << endl;

  Vector2cd psi = initialState.normalized();
  for (double t = tmin; t < tmax; t += dt) {
    Vector2cd phi = Phi(psi, t), phin = phi.normalized();
    out << real((projector(psi)*sigma_x_pm).trace()) << " " << real((projector(phin)*sigma_x_pm).trace()) << " " << real((projector(psi)*sigma_z_pm).trace()) << " " << real((projector(phin)*sigma_z_pm).trace()) << " " << phi.norm() << endl;
    psi += -I*dt*(H(t) - .5*I*Gamma(t))*psi - .5*dt*phi;
    psi.normalize();
  }
}