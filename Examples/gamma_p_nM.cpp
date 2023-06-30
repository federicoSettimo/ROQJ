// Non-Markovianity from gamma_p < 0
// Writing everyting in the |+,-> basis, hoping to avoid |0>
#include "../roqj.h"

using namespace std;
using namespace Eigen;

Matrix2cd sigma_p_pm {{1,-1},{1,-1}}, sigma_m_pm {{1,1},{-1,-1}}, sigma_z_pm {{0,1},{1,0}}, sigma_x_pm {{1,0},{0,-1}}, sigma_y_pm {{0,I},{-I,0}};

// Note: it must be gp + gm + gz >= 0 at all times
double kappa = .35;
/*double gamma_p (double t) {return exp(-t)*(kappa + (1.-kappa)*exp(-.25*t)*cos(2.*t));}
double gamma_m (double t) {return exp(-2.*t);}
double gamma_z (double t) {return exp(-3.*t);}*/
// Second choice of functions: gp big for t small => to northern hemisphere
double sig(double t) {return 1./(1.+exp(-15.*t));}
double gamma_p (double t) {return sig(1.-t) + sig(t-1.)*exp(-t)*(kappa + (1.-kappa)*exp(-.25*t)*cos(2.*t));}
double gamma_m (double t) {return exp(-7.*t);}
double gamma_z (double t) {return exp(-.25*t);}


MatrixXcd H (double t) {
  return MatrixXcd::Zero(2,2);
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return gamma_p(t)*sigma_p_pm*rho*sigma_m_pm + gamma_m(t)*sigma_m_pm*rho*sigma_p_pm + gamma_z(t)*sigma_z_pm*rho*sigma_z_pm;
}

MatrixXcd Gamma (double t) {
  return gamma_p(t)*sigma_m_pm*sigma_p_pm + gamma_m(t)*sigma_p_pm*sigma_m_pm + gamma_z(t)*id;
}

MatrixXcd C (const VectorXcd &psi, double t) {
  Vector2cd cpsi = psi;
  if (real(psi(0)) < 0.) cpsi = -cpsi;
  double a = real(cpsi(1)), cp, cm; // a = <psi|->
  Matrix2cd CC = MatrixXcd::Zero(2,2);
  double a2 = norm(a), gm = gamma_m(t), gp = gamma_p(t), gz = gamma_z(t);

  if (a2 >= .99) {// |->
    CC(0,1) = 2.*(gm-gp);
    CC(1,1) = gz;
  }
  else if (a2 <= .01) { // |+>
    CC(0,0) = gz;
    CC(1,0) = 2.*(gm-gp);
  }
  else { // ok unless |0> (or close), but hoping to avoid it
    double cub, clb;
    cub = 2.*(gm-gp)*(a+1./a)/sqrt(1.-a2) + (5.*gp+5.*gm-2.*gz)/(1.-a2) + (3.*gz-4.*gp-4.*gm)*a2/(1.-a2);
    clb = 2.*(gm-gp)*sqrt(1.-a2)/a + (gm+gp)/a2 + gz*(1.-a2)/a2;
    //cm = a2*clb + (1.-a2)*cub;
    cm = cub;
    cp = -cm + 4.*(gp+gm) - 2.*gz + 2.*(gm-gp)/(a*sqrt(1.-a2));

    CC(0,0) = cp;
    CC(1,1) = cm;
  }
  return CC;
}

/*MatrixXcd C (const VectorXcd &psi, double t) {
  Vector2cd cpsi = psi, psi_perp;
  psi_perp(0) = cpsi(1);
  psi_perp(1) = -cpsi(0);
  if (real(psi(0)) < 0.) cpsi = -cpsi;
  double a = real(cpsi(1)), cp, cm; // a = <psi|->
  double a2 = norm(a), a4 = a2*a2, sq = sqrt(1.-a2), gm = gamma_m(t), gp = gamma_p(t), gz = gamma_z(t);
  double c,d;

  if (abs(a) <= .5/sqrt(2.)) { // around 0
    c = -((1.-6.*a2+4.*a4-4.*a*a2*sq)*gm + (1.-6.*a2+4.*a4+4.*a*a2*sq)*gp + (3.*a2-4.*a4)*gz)/(1.-a2);
    d = 2.*(c*a*sq - (1.+2.*a*sq)*gm + (1.-2.*a*sq)*gp + a*sq*gz)/(1.-2.*a2);
  }
  else if (abs(a) >= .5*(1.+1./sqrt(2.))) { // around +/-1
    c = -((-1.-2.*a2+4.*a4-4.*a*sq+4.*a*a2*sq)*gm + (-1.-2.*a2+4.*a4+4.*a*sq-4.*a*a2*sq)*gp + (-1.+5.*a2-4.*a4))/a2;
    d = 2.*(c*a*sq - (1.+2.*a*sq)*gm + (1.-2.*a*sq)*gp + a*sq*gz)/(1.-2.*a2);
  }
  else {
    double dmub, dmlb;
    dmub = 2.*(1.+2.*a2)*(gp-gm) + a*(2.*gz-6.*gp-6.*gm)/sq + 4.*a*a2*(gz-gp-gm)/sq;
    dmlb = 2.*sq*(gm+gp+gz)/a + 2.*(3.-2.*a2)*(gm-gp) + 4.*a*sq*(gp+gm-gz);
    d = .5*(dmub+dmlb);
    c = .5*((1.-2.*a2)*d + 2.*(1.+2.*a*sq)*gm + 2.*(-1.+2.*a*sq)*gp - 2.*a*sq*gz)/(a*sq);
  }

  return c*cpsi*cpsi.transpose() + d*psi_perp*cpsi.transpose();
}*/

double observable (const MatrixXcd &rho) {return real((rho*sigma_z_pm).trace());}

int main () {
  double tmin = 0., tmax = 3., dt = 0.01, threshold = 1e-10;
  int N_ensemble = 1000, Ncopies = 1, dimH = 2, Ntraj = 10;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj, true, threshold);
  Vector2cd initialState;
  initialState << 0.2, 0.6;
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}