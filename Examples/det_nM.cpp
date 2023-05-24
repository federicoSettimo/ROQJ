/*
  Comparison between the total quantum nM and the nM arising only from the det. evolution
  Not really working at the moment
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

MatrixXcd C (const Vector2cd &psi, double t, double kappa, double theta) {
  double mu, c3, a2 = norm(psi(1)), eps = -.5*sqrt(gamma_m(t)*gamma_p(t)) - gamma_z(t,kappa), mu_ub, mu_lb;
  eps = eps > 0. ? eps : 0.;
  if (theta < 1.3)
    mu = a2 == 0. ? sqrt(gamma_p(t)*gamma_m(t)) + 2.*eps : 2.*gamma_z(t,kappa) + gamma_m(t)*(1-a2)/a2;
  else mu = a2 == 1. ? -sqrt(gamma_m(t)*gamma_p(t)) - 2.*eps : ((1.-a2)*sqrt(gamma_p(t)*gamma_m(t)) - a2*gamma_p(t))/(1.-a2) + 2.*eps;

  //mu = a2 == 0. ? sqrt(gamma_p(t)*gamma_m(t)) + 2.*eps : 2.*gamma_z(t,kappa) + gamma_m(t)*(1-a2)/a2;

  c3 = gamma_z(t,kappa) - mu;
  return 2.*mu*sigma_p*sigma_m + c3*id;
}

Matrix2cd projector (const Vector2cd &psi) {return psi*psi.adjoint();}

Matrix2cd comm (const Matrix2cd &A, const Matrix2cd &B) {return A*B - B*A;}
Matrix2cd anticomm (const Matrix2cd &A, const Matrix2cd &B) {return A*B + B*A;}

double get_non_Markovianity (double kappa);

double get_non_Markovianity_det (double kapppa);

double TD (const Matrix2cd &rho, const Matrix2cd &sigma);

int main () {
  double dkappa = 0.05;

  ofstream out, np;
  out.open("det_nM.txt");
  np.open("det_nM_np.txt");

  int npoints = 0.;
  for (double kappa = 1.; kappa < 8.25; kappa += dkappa) {
    cout << "kappa = " << kappa << endl;
    out << kappa << " " << get_non_Markovianity(kappa) << " " << get_non_Markovianity_det(kappa) << endl;
    npoints++;
  }
  np << npoints;

  return 0;
}

double dtheta = 0.05, dt = 0.01, tmax = 10.;

double get_non_Markovianity (double kappa) {
  double Nmax = 0.;
  for (double theta = 0.; theta < .5*M_PI + dtheta; theta += dtheta) {
    if (theta > .5*M_PI) theta = .5*M_PI;
    Vector2cd psi0 = cos(theta)*ground_state + sin(theta)*excited_state, psi_perp = cos(M_PI-theta)*ground_state + sin(M_PI-theta)*excited_state;
    Matrix2cd rho1 = projector(psi0), rho2 = projector(psi_perp);

    double oldd = 1., newd = 1., N = 0.;
    for (double t = 0.; t < tmax; t += dt) {
      oldd = newd;
      rho1 += (J(rho1,t,kappa) - .5*anticomm(Gamma(t,kappa),rho1))*dt;
      rho2 += (J(rho2,t,kappa) - .5*anticomm(Gamma(t,kappa),rho2))*dt;

      newd = TD(rho1,rho2);
      if (newd > oldd) N += (newd-oldd);
    }

    if (N > Nmax) Nmax = N;
  }
  return Nmax;
}

/*
double get_non_Markovianity_det (double kappa) {
  double Nmax = 0.;
  for (double theta = 0.; theta < .5*M_PI + dtheta; theta += dtheta) {
    if (theta > .5*M_PI) theta = .5*M_PI;
    Vector2cd psi1 = cos(theta)*ground_state + sin(theta)*excited_state, psi2 = cos(M_PI-theta)*ground_state + sin(M_PI-theta)*excited_state;

    double oldd = 1., newd = 1., N = 0., pdet1 = 1., pdet2 = 1.;
    for (double t = 0.; t < tmax; t += dt) {
      oldd = newd;
      
      psi1 -= .5*(Gamma(t, kappa) + C(psi1, t, kappa, theta))*psi1*dt;
      psi2 -= .5*(Gamma(t, kappa) + C(psi2, t, kappa, M_PI-theta))*psi2*dt;
      psi1.normalize();
      psi2.normalize();

      Matrix2cd rho1 = projector(psi1), rho2 = projector(psi2), R = J(rho1, t, kappa) + .5*(C(psi1,t,kappa,theta)*rho1 + rho1*C(psi1,t,kappa,theta).adjoint());
      pdet1 -= dt*real(R.trace());
      if (real(R.trace()) < 0.) break;
      R = J(rho2, t, kappa) + .5*(C(psi2,t,kappa,M_PI-theta)*rho2 + rho2*C(psi2,t,kappa,M_PI-theta).adjoint());
      pdet2 -= dt*real(R.trace());
      if (real(R.trace()) < 0.) break;

      newd = TD(pdet1*rho1,pdet2*rho2);
      if (newd > oldd) N += (newd-oldd);
    }

    if (N > Nmax) Nmax = N;
  }
  return Nmax;
}
*/

/*
double get_non_Markovianity_det (double kappa) {
  double Nmax = 0.;

  for (double theta1 = 0.; theta1 < M_PI + dtheta; theta1 += dtheta) {
    bool ok_theta1 = true;

    // Cycle also on theta2. If theta1 gives some <0 rates, just skip this theta1
    for (double theta2 = 0.; theta2 < .5*M_PI && ok_theta1; theta2 += dtheta) {
      bool ok_theta2 = true;
      Vector2cd psi1 = cos(theta1)*ground_state + sin(theta1)*excited_state, psi2 = cos(theta2)*ground_state + sin(theta2)*excited_state;
      psi1.normalize();
      psi2.normalize();

      double oldd = 1., newd = 1., N = 0., pdet1 = 1., pdet2 = 1.;
      for (double t = 0.; t < tmax && ok_theta1 && ok_theta2; t += dt) {
        // Cycle on t, but only for the initial states giving a positive unravelling at all times
        oldd = newd;

        psi1 -= .5*(Gamma(t, kappa) + C(psi1, t, kappa, theta1))*psi1*dt;
        psi2 -= .5*(Gamma(t, kappa) + C(psi2, t, kappa, theta2))*psi2*dt;
        psi1.normalize();
        psi2.normalize();

        Matrix2cd rho1 = projector(psi1), rho2 = projector(psi2), R = J(rho1, t, kappa) + .5*(C(psi1,t,kappa,theta1)*rho1 + rho1*C(psi1,t,kappa,theta1).adjoint());
        pdet1 -= dt*real(R.trace())*pdet1;
        if (real(R.trace()) < 0.) ok_theta1 = false;
        R = J(rho2, t, kappa) + .5*(C(psi2,t,kappa,theta2)*rho2 + rho2*C(psi2,t,kappa,theta2).adjoint());
        pdet2 -= dt*real(R.trace())*pdet2;
        if (real(R.trace()) < 0.) ok_theta2 = false;

        newd = TD(pdet1*rho1,pdet2*rho2);
        if (newd > oldd) N += (newd-oldd);
      }
      if (N > Nmax)
        Nmax = N;
    }
  }
  return Nmax;
}
*/

double get_non_Markovianity_det (double kappa) {
  double Nmax = 0.;

  for (double theta = 0.; theta < .5*M_PI; theta += dtheta) {
    bool theta_ok = true;
    Vector2cd psi = (cos(theta)*ground_state + sin(theta)*excited_state).normalized();

    double oldd = 1., newd = 1., N = 0., pdet = 1.;
    for (double t = 0.; t < tmax && theta_ok; t += dt) {
      psi -= .5*(Gamma(t, kappa) + C(psi, t, kappa, theta))*psi*dt;
      psi.normalize();

      Matrix2cd rho = projector(psi), R = J(rho, t, kappa) + .5*(C(psi,t,kappa,theta)*rho + rho*C(psi,t,kappa,theta).adjoint());

      double dpdet = real(R.trace());
      if (dpdet < -1e-10) theta_ok = false;
      pdet -= dt*dpdet*pdet;

      oldd = newd;
      newd = 2.*real(projector(psi)(0,1))*pdet; // Revival in TD = revival in coherence
      if (newd > oldd && theta_ok) {
        //cout << t << " " << kappa << " " << theta << ", " << newd << " " << oldd << " " << newd - oldd << " " << pdet << endl;
        N += (newd - oldd);
      }
    }
    if (N > Nmax)
      Nmax = N;
  }
  return Nmax;
}

double TD (const Matrix2cd &rho, const Matrix2cd &sigma) {
  Matrix2cd A = rho - sigma;
  Matrix2cd B = A.transpose()*A;
  return .5*real((B.sqrt()).trace());
}