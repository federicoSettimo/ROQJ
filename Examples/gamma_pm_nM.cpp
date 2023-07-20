// Non-Markovianity from gamma_p < 0
// Writing everyting in the |+,-> basis, hoping to avoid |0>
#include "../roqj.h"

using namespace std;
using namespace Eigen;

Matrix2cd sigma_p_pm {{1,-1},{1,-1}}, sigma_m_pm {{1,1},{-1,-1}}, sigma_z_pm {{0,1},{1,0}}, sigma_x_pm {{1,0},{0,-1}}, sigma_y_pm {{0,I},{-I,0}};

// Note: it must be gp + gm + gz >= 0 at all times
double kappa = .35;
double gamma_p (double t) {return exp(-.1*t)*(kappa + (1.-kappa)*exp(-.25*t)*cos(2.*t));}
double gamma_m (double t) {return gamma_p(t);}
double gamma_z (double t) {return 1.;}

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
  if (real(psi(0)) < 0.) cpsi -= psi;
  double a = real(cpsi(1)), cp, cm; // a = <psi|->
  Matrix2cd CC = MatrixXcd::Zero(2,2);
  double a2 = norm(a), g = gamma_m(t), gz = gamma_z(t);

  double cub, clb, c;
  //cub = 2.*g*(1.+4.*a2)/a2 + gz*(1.-2.*a2)/a2;
  //clb = -(2.*g + gz)/(1.-a2);

  clb = (2.*g + a2*gz)/(-1 + a2);
  cub = (2.*g + 8.*a2*g + gz - 3.*a2*gz)/a2;
  if (a2 <= .01) c = clb;
  else if (a2 >= .99) c = cub;
  else  if (a > 1./sqrt(2.)) c = cub;
  else c = clb;

  CC(0,0) = c;
  CC(1,1) = 8.*g - 2.*gz - c;

  return CC;
}

double observable (const MatrixXcd &rho) {return real((rho*sigma_z_pm).trace());}

int main () {
  double tmin = 0., tmax = 2., dt = 0.0005, threshold = 1e-3;
  int N_ensemble = 500, Ncopies = 3, dimH = 2, Ntraj = 15;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj, true, threshold);
  Vector2cd initialState;
  initialState << 0.2, -0.3;
  //initialState << 1,-1;
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  // Plotting the functions
  ofstream out;
  out.open("functions.txt");
  for (double t = tmin; t < tmax; t += dt)
    out << gamma_p(t) << " " << gamma_m(t) << " " << gamma_z(t) << endl;

  return 0;
}