// Phase covariant violating also P divisibility - no reverse jumps
// Dynamics from Tetteinen - Quantum speed limit and divisibility of the dynamical map
// Tracking with only 1 bit: jumps |psi(t)> -> |0> and |0> -> |psi(t+dt)>
#include "../roqj_state.h"

using namespace std;
using namespace Eigen;

Matrix2cd sigma_p_pm {{1,-1},{1,-1}}, sigma_m_pm {{1,1},{-1,-1}}, sigma_z_pm {{0,1},{1,0}}, sigma_x_pm {{1,0},{0,-1}}, sigma_y_pm {{0,I},{-I,0}};
Vector2cd plus_state_pm {{1,0}}, minus_state_pm {{0,1}};

double tmin = 0., tmax = 3., dt = 0.0005, threshold = 1e-3;
vector<Vector2cd> psi_t((int)((tmax-tmin)/dt+1));

double kappa = .5;
double gamma_p (double t) {return exp(-.1*t)*(kappa + (1.-kappa)*exp(-.25*t)*cos(2.*t));}
double gamma_m (double t) {return gamma_p(t);}
double gamma_z (double t) {return .5;}

MatrixXcd H (double t) {
  return MatrixXcd::Zero(2,2);
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return gamma_p(t)*sigma_p_pm*rho*sigma_m_pm + gamma_m(t)*sigma_m_pm*rho*sigma_p_pm + gamma_z(t)*sigma_z_pm*rho*sigma_z_pm;
}

MatrixXcd Gamma (double t) {
  return gamma_p(t)*sigma_m_pm*sigma_p_pm + gamma_m(t)*sigma_p_pm*sigma_m_pm + gamma_z(t)*id;
}

VectorXcd Phi (const VectorXcd &psi, double t) {
  double a = real(psi(1)), a2 = a*a, sq = sqrt(1.-a2), g = gamma_p(t), gz = gamma_z(t);

  // From |psi>
  if (a2 < 1.)
    return (-2.*g - a2*gz)/sq*plus_state_pm - a*(2.*(4.*a2-5.)*g + (2.-3.*a2)*gz)/(sq*sq)*minus_state_pm;

  // From |0>
  Vector2cd target = psi_t[(int)(t/dt)+1]; // Target = state at the next timestep
  double b = real(target(1)), b2 = b*b, sqb = sqrt(1.-b*b);
  return 2.*b*(2.*g+gz)/sqb*plus_state_pm - (2.*(1.-2.*b2)*g - b2*gz)/(sqb*sqb)*minus_state_pm;
}

double observable (const MatrixXcd &rho) {return real((rho*sigma_z_pm).trace());}

int main () {
  int N_ensemble = 10000, Ncopies = 3, dimH = 2, Ntraj = 10;
  bool printTraj = true;

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj, true, threshold);
  Vector2cd initialState;
  initialState << 0.7, 0.6;
  jump.set_initial_state(initialState);

  // First thing: generating the deterministic evolution
  int i = 1;
  psi_t[0] = initialState.normalized();
  for (double t = tmin; t < tmax; t += dt) {
    Matrix2cd K = H(t) - .5*I*Gamma(t);
    psi_t[i] = psi_t[i-1] - I*dt*K*psi_t[i-1] - .5*dt*Phi(psi_t[i-1],t);
    psi_t[i].normalize();
    i++;
  }

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}