// Phase covariant violating P divisibility - without reverse jumps
// Jumps to only 1 post-jump state (fixed)
#include "../roqj.h"

using namespace std;
using namespace Eigen;

double gamma_p (double t) {return exp(-.5*t);}
double gamma_m (double t) {return exp(-.25*t);}
double gamma_z (double t) {return 1.3*exp(-3.*t/8.)*cos(2.*t)*.5;}

double alpha0 = .1;

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

MatrixXcd C (const VectorXcd &psi, double t) {
  double theta = acos(alpha0), ag = abs(psi.dot(ground_state)), ae = sqrt(1.-ag*ag), gz = gamma_z(t), gp = gamma_p(t), gm = gamma_m(t);

  double eps = -.5*sqrt(gm*gp) - gz;
  eps = eps > 0. ? eps : 0.;

  if (abs(ag) < 0.01) {
    double mu = sqrt(gamma_p(t)*gamma_m(t)) + 2.*eps;
    double c3 = gamma_z(t) - mu;
    return 2.*mu*sigma_p*sigma_m + c3*id;
  }

  if (abs(1.-ag) < 0.01) {
    double mu = 2.*gz;
    double c3 = gamma_z(t) - mu;
    return 2.*mu*sigma_p*sigma_m + c3*id;
  }

  MatrixXcd CC = MatrixXcd::Zero(2,2);
  if (abs(ag-alpha0) < .01) theta = M_PI-theta; // If close to the post-jump state, jumps to its orthogonal

  double c1 = (-(ae*pow(ag,2)*(gp - 3.*gz)) + pow(ae,3)*(gm - gz) - ae*(pow(ae,2)*(gm + gz) + pow(ag,2)*(gp + 3.*gz))*cos(2.*theta) + 
     2*ag*(pow(ag,2)*gp + pow(ae,2)*gz)*sin(2*theta))/(4.*ae*pow(ae*cos(theta) - ag*sin(theta),2));
  double c2 = -0.25*(pow(ae,2)*ag*(gm - 3.*gz) + pow(ag,3)*(-gp + gz) - ag*(pow(ag,2)*(gp + gz) + pow(ae,2)*(gm + 3.*gz))*cos(2.*theta) - 
      2*ae*(pow(ae,2)*gm + pow(ag,2)*gz)*sin(2.*theta))/(ag*pow(ae*cos(theta) - ag*sin(theta),2));

  CC(0,0) = c1;
  CC(1,1) = c2;
  return CC;
}

double observable (const MatrixXcd &rho) {return abs(rho(0,1));}
//double observable (const MatrixXcd &rho) {return real((rho*sigma_z).trace());}

int main () {
  double tmin = 0., tmax = 3., dt = 0.01, threshold = 1e-3;
  int N_ensemble = 10000, Ncopies = 1, dimH = 2, Ntraj = 10;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj, true, threshold);
  Vector2cd initialState;
  initialState << 0.1, 0.9;
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}