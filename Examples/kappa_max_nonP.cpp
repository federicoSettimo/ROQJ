/*
  Non-P divisible dynamics [Tetteinen - Quantum speed limit and divisibility of the dynamical map],
  checking the initial states for which a given value of kappa does not give <0 eigenvalues of R
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

MatrixXcd C (const Vector2cd &psi, double t, double kappa) {
  double mu, c3, a2 = norm(psi(1)), eps = -.5*sqrt(gamma_m(t)*gamma_p(t)) - gamma_z(t,kappa);
  eps = eps > 0. ? eps : 0.;
  mu = a2 == 0. ? sqrt(gamma_p(t)*gamma_m(t)) + 2.*eps : 2.*gamma_z(t,kappa) + gamma_m(t)*(1-a2)/a2;
  c3 = gamma_z(t,kappa) - mu;
  return 2.*mu*sigma_p*sigma_m + c3*id;
}

bool isPdiv (double t, double kappa) {return sqrt(gamma_p(t)*gamma_m(t)) + 2.*gamma_z(t,kappa) >= -1e-5;}

Matrix2cd projector (const Vector2cd &psi) {return psi*psi.adjoint();}

int main () {
  double dkappa = .001, dtheta = .0001*M_PI, dt = .01, threshold = 1e-13;
  int npoints = 0; // Number of points (kappa,theta) printed

  ofstream out, npoints_out;
  out.open("kappa_theta.txt");
  npoints_out.open("kappa_theta_npoints.txt");

  // K = 8.25 is approximately the threshold after which the map is no more CP for all times
  for (double kappa = 1.; kappa <= 8.25; kappa += dkappa) {
    for (double theta = 0.; theta < .5*M_PI + dtheta; theta += dtheta) { // theta = 0: |g>; theta = pi/2: |e>
      Vector2cd initialState = cos(theta)*ground_state + sin(theta)*excited_state, psi;
      initialState.normalize();
      if (theta >= .5*M_PI) initialState = excited_state;
      psi = initialState;

      // Checks whether nM has begun or ended to stop the time-evol. theta_ok = true iff that theta has >0 eigs
      bool finished_nM = false, started_nM = false, theta_ok = true;
      for (double t = 0.; t <= 100. && !finished_nM && theta_ok; t += dt) {
        Matrix2cd K = .5*(C(psi,t,kappa).imag() - I*(Gamma(t,kappa) + C(psi,t,kappa).real())), R;
        psi -= K*psi*I*dt;
        psi.normalize();

        // If P divisible, do nothing (it is already ok)
        if (!isPdiv(t,kappa)) {
          if(!started_nM)
            started_nM = true;
          R = J(projector(psi),t,kappa) + 0.5*(C(psi, t,kappa)*projector(psi) + projector(psi)*(C(psi, t,kappa).adjoint()));
          ComplexEigenSolver<Matrix2cd> eigs;
          eigs.compute(R);
          Vector2cd eigval = eigs.eigenvalues();
          if (real(eigval(0)) < -threshold || real(eigval(1)) < -threshold)
            theta_ok = false;
        }
        else if (started_nM)
          finished_nM = true;
      }

      if (theta_ok) {
        out << kappa << " " << theta << endl;
        npoints++;
      }
    }
  }
  npoints_out << npoints << endl;
  return 0;
}