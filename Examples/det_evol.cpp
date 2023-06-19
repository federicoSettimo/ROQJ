/*
  Plotting the deterministic evolution of the states and the evolution of its norm
*/

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
#include <unsupported/Eigen/MatrixFunctions>
#include <cstring>

using namespace std;
using namespace Eigen;

complex<double> I(0,1), one(1,0), zero(0,0);
Eigen::Matrix2cd sigma_x {{0,1},{1,0}}, sigma_y {{0,-I},{I,0}}, sigma_z {{1,0},{0,-1}}, sigma_p {{0,1},{0,0}}, sigma_m {{0,0},{1,0}}, id {{1,0},{0,1}};
Eigen::Vector2cd ground_state {{0.,1.}}, excited_state {{1.,0.}}, plus_state {{1./sqrt(2.),1./sqrt(2.)}}, minus_state {{1./sqrt(2.),-1./sqrt(2.)}};

double gamma_p (double t) {return exp(-.5*t);}
double gamma_m (double t) {return exp(-.25*t);}
double gamma_z (double t, double kappa) {return .5*kappa*exp(-3.*t/8.)*cos(2.*t);}

MatrixXcd J (const MatrixXcd &rho, double t, double kappa) {
  return gamma_p(t)*sigma_p*rho*sigma_m + gamma_m(t)*sigma_m*rho*sigma_p + gamma_z(t,kappa)*sigma_z*rho*sigma_z;
}

MatrixXcd Gamma (double t, double kappa) {
  return gamma_p(t)*sigma_m*sigma_p + gamma_m(t)*sigma_p*sigma_m + gamma_z(t,kappa)*id;
}

double alpha0 = .2;

MatrixXcd C (const Vector2cd &psi, double t, double kappa, double theta) {
  MatrixXcd CC = MatrixXcd::Zero(2,2);
  double ag = abs(psi.dot(ground_state)), ae = sqrt(1.-ag*ag), gz = gamma_z(t,kappa), gp = gamma_p(t), gm = gamma_m(t);

  if (ag == alpha0) {
    double c1 = (gm - 7.*gp + 4*gz - 8.*(gp - 2.*gz)*cos(2.*theta) - (gm + gp - 4.*gz)*cos(4*theta))/16.;
    double c2 = (-7.*gm + gp + 4.*gz + 8.*(gm - 2.*gz)*cos(2.*theta) - (gm + gp - 4.*gz)*cos(4*theta))/16.;
    CC(0,0) = c1;
    CC(1,1) = c2;
    return CC;
  }

  double c1 = (-(ae*pow(ag,2)*(gp - 3.*gz)) + pow(ae,3)*(gm - gz) - ae*(pow(ae,2)*(gm + gz) + pow(ag,2)*(gp + 3.*gz))*cos(2.*theta) + 
     2*ag*(pow(ag,2)*gp + pow(ae,2)*gz)*sin(2*theta))/(4.*ae*pow(ae*cos(theta) - ag*sin(theta),2));
  double c2 = -0.25*(pow(ae,2)*ag*(gm - 3.*gz) + pow(ag,3)*(-gp + gz) - ag*(pow(ag,2)*(gp + gz) + pow(ae,2)*(gm + 3.*gz))*cos(2.*theta) - 
      2*ae*(pow(ae,2)*gm + pow(ag,2)*gz)*sin(2.*theta))/(ag*pow(ae*cos(theta) - ag*sin(theta),2));

  
  CC(0,0) = c1;
  CC(1,1) = c2;
  return CC;
}

Matrix2cd projector (const Vector2cd &psi) {return psi*psi.adjoint();}

int main (int argc, char *argv[]) {
  ofstream out;
  out.open("det_evol.txt");

  double kappa = 1.3, theta = acos(alpha0); // Theta = angle of the fixed state
  double tmax = 4., dt = 0.01, dalpha = 0.005;
  int Nsteps = (int)(tmax/dt)+1, Nstates = (int)(1./dalpha)-1;
  out << Nsteps << endl << Nstates << endl << alpha0 << endl;

  // For 0,1 not ok: modify later
  for (double alpha = dalpha; alpha < 1.; alpha += dalpha) {
    Vector2cd psi = (alpha*ground_state + sqrt(1.-alpha*alpha)*excited_state).normalized();
    for (double t = 0.; t < tmax; t += dt) {
      out << t << " " << abs(real(psi.dot(ground_state))) << endl;
      psi -= .5*(Gamma(t, kappa) + C(psi.normalized(), t, kappa, theta))*psi*dt;
      psi.normalize();
    }
  }

  // Printing the not-ok-alpha region
  ofstream alpha_no;
  alpha_no.open("alpha_no.txt");
  double t_noP_min = tmax, t_noP_max;
  double alpha_min = 1., alpha_max = 0.;
  for (double t = 0.; t < tmax; t += dt) {
    double gp = gamma_p(t), gm = gamma_m(t), gz = gamma_z(t,kappa), Pdiv = gz + .5*sqrt(gm*gp);
    if (Pdiv < 0.) {
      double alpha =  1./sqrt(1. + sqrt(gm/gp));
      if (alpha < alpha_min)
        alpha_min = alpha;
      if (alpha > alpha_max)
        alpha_max = alpha;
      if (t < t_noP_min)
        t_noP_min = t;
      t_noP_max = t;
    }
  }
  alpha_no << alpha_min << endl << alpha_max << endl << t_noP_min << endl << t_noP_max;

  // Now det evolution but starting from a jump at some time tjump
  ofstream out_jump;
  out_jump.open("det_evol_jump.txt");
  Vector2cd psi0 = (alpha0*ground_state + sqrt(1.-alpha0*alpha0)*excited_state).normalized(), psi;
  for (double tjump = 0.; tjump < tmax; tjump += dt) {
    psi = psi0;
    for (double t = 0.; t < tmax; t += dt) {
      out_jump << t << " " << abs(real(psi.dot(ground_state))) << endl;
      if (t > tjump) {
        psi -= .5*(Gamma(t, kappa) + C(psi.normalized(), t, kappa, theta))*psi*dt;
        psi.normalize();
      }
    }
  }
  // Same but with the orthogonal (2 jumps)
  psi0 = (sqrt(1.-alpha0*alpha0)*ground_state - alpha0*excited_state).normalized();
  for (double tjump = 0.; tjump < tmax; tjump += dt) {
    psi = psi0;
    for (double t = 0.; t < tmax; t += dt) {
      out_jump << t << " " << abs(real(psi.dot(ground_state))) << endl;
      if (t > tjump) {
        psi -= .5*(Gamma(t, kappa) + C(psi.normalized(), t, kappa, theta))*psi*dt;
        psi.normalize();
      }
    }
  }

  return 0;
}
