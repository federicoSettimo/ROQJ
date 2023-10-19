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
Matrix2cd sigma_x {{0,1},{1,0}}, sigma_y {{0,-I},{I,0}}, sigma_z {{1,0},{0,-1}}, sigma_p {{0,1},{0,0}}, sigma_m {{0,0},{1,0}}, id {{1,0},{0,1}};

Vector2cd ground_state {{0.,1.}}, excited_state {{1.,0.}}, plus_state {{1./sqrt(2.),1./sqrt(2.)}}, minus_state {{1./sqrt(2.),-1./sqrt(2.)}};

Matrix2cd projector (const Vector2cd &psi) {return psi*psi.adjoint();}

double gamma_m (double t) {return 1.;}
double gamma_p (double t) {return gamma_m(t);}
double beta (double t) {return 2.;}

Matrix2cd H (double t) {
  return beta(t)*sigma_y;
}

Matrix2cd J (const Matrix2cd &rho, double t) {
  return gamma_m(t)*sigma_m*rho*sigma_p + gamma_p(t)*sigma_p*rho*sigma_m;
}

Matrix2cd Gamma (double t) {
  return gamma_m(t)*sigma_p*sigma_m + gamma_p(t)*sigma_m*sigma_p;
}

double observable (const Matrix2cd &rho) {
  return real((rho*sigma_z).trace());
}

double coherence (const Matrix2cd &rho) {
  return real((rho*sigma_x).trace());
}

Matrix2cd comm (const Matrix2cd &A, const Matrix2cd &B) {return A*B - B*A;}

Matrix2cd anticomm (const Matrix2cd &A, const Matrix2cd &B) {return A*B + B*A;}

Vector3d BlochVector (const Matrix2cd &A) {
  Vector3d r;
  r << real((A*sigma_x).trace()), real((A*sigma_y).trace()), real((A*sigma_z).trace());
  return r;
}

double TD (const Matrix2cd &A, const Matrix2cd &B) {return .5*(BlochVector(A) - BlochVector(B)).norm();}

Vector2cd Phi (double a, double t) {
  double g = gamma_p(t);
  return (2*pow(a,2)*(1 - 2*pow(a,2))*sqrt(1 - pow(a,2))*g)*excited_state + (-2*a*(1 - 3*pow(a,2) + 2*pow(a,4))*g)*ground_state;
}

Vector2cd PhiPerp (double a, double t) {
  double g = gamma_p(t);
  return (2*a*(1 - 3*pow(a,2) + 2*pow(a,4))*g)*excited_state + (2*pow(a,2)*(1 - 2*pow(a,2))*sqrt(1 - pow(a,2))*g)*ground_state;
}

Vector2cd orthogonal (const Vector2cd &psi) {
  double a = real(psi(1));
  return -a*excited_state + sqrt(1.-a*a)*ground_state;
}

double tmin = 0., tmax = 5., dt = .0001, threshold = 1e-3;
int Ncopies = 10, Nensemble = 100, Npsi = Nensemble, NPerp = 0;
int Npsi_old, NPerp_old;

void jump (double p,int &N, int &N1) {
  double z = (double)rand()/(double)RAND_MAX;
  if (z < p) {N--; N1++;}
}

Matrix2cd L (const Matrix2cd &rho, double t) {return -I*comm(H(t),rho) + J(rho,t) - .5*anticomm(Gamma(t),rho);}

// Runge-Kutta 4
Matrix2cd RK4 (const Matrix2cd &rho, double t) {
  Matrix2cd k1 = dt * L(rho, t), k2 = dt * L(rho + .5*k1, t + .5*dt), k3 = dt * L(rho + .5*k2, t + .5*dt), k4 = dt * L(rho + k3, t + dt);
  return (1./6.) * (k1 + 2.*k2 + 2.*k3 + k4);
}

int main () {
  srand(time(NULL));
  ofstream out;
  out.open("driven.txt");
  out << (int)(tmax/dt) << endl << tmax << endl << dt << endl;

  Vector2cd initialState = .7*excited_state + .8*ground_state, psi, psiPerp;
  initialState.normalize();
  psi = initialState;
  psiPerp = orthogonal(psi);

  int Nsteps = (int)((tmax)/dt) + 1;
  MatrixXd observables(Ncopies, Nsteps), coherences(Ncopies, Nsteps), errors(Ncopies, Nsteps);
  Matrix2cd K;

  cout << "Using " << Ncopies << " copies, each with " << Nensemble << " states. dt = " << dt << endl;
  psi = initialState.normalized();
  Matrix2cd exact = projector(psi), rho = projector(psi);
  for (int Ncopy = 0; Ncopy < Ncopies; ++Ncopy) {
    cout << "Running copy " << Ncopy+1 << "/" << Ncopies << "...\n";

    psi = initialState.normalized(); exact = projector(psi); psiPerp = orthogonal(psi);
    int Nstep = 0;
    Npsi = Nensemble; NPerp = 0;

    for (double t = 0.; t < tmax; t += dt) {
      // Old populations
      Npsi_old = Npsi; NPerp_old = NPerp;
      // Rates and parameters
      double a = real(psi(1)), gp = gamma_p(t), gm = gamma_m(t), a2 = a*a, b = beta(t);
      // Fractions of states
      double fpsi = (double)Npsi/(double)Nensemble, fPerp = (double)NPerp/(double)Nensemble;      
      rho = fpsi*projector(psi) + fPerp*projector(psiPerp);
      if (Ncopy == 0) {
        out << observable(exact) << " " << observable(projector(psi)) << endl;
        out << coherence(exact) << " " << coherence(projector(psi)) << endl;
        out << fpsi << " " << fPerp << endl;
        out << gp << " " << gm << " " << b << endl;
      }
      observables(Ncopy, Nstep) = observable(rho);
      coherences(Ncopy, Nstep) = coherence(rho);
      errors(Ncopy, Nstep) = TD(rho, exact);
      
      // Updating exact
      exact += RK4(exact,t);

      // Updating psi
      Vector2cd phi = Phi(a, t), phiPerp = PhiPerp(a, t);
      Matrix2cd R = J(projector(psi), t) + .5*(psi*phi.adjoint() + phi*psi.adjoint()), RPerp = J(projector(psiPerp), t) + .5*(psiPerp*phiPerp.adjoint() + phiPerp*psiPerp.adjoint());
      Matrix2cd K = H(t) - .5*I*Gamma(t);
      psi -= dt*(I*K*psi + .5*phi); psi.normalize();
      psiPerp -= dt*(I*K*psiPerp + .5*phiPerp); psiPerp.normalize();
      
      double psi_psiPerp = 0., psiPerp_psi = 0., lpsi_psiPerp = 0., lpsiPerp_psi = 0.;
      ComplexEigenSolver<Matrix2cd> eigs;
      eigs.compute(R);
      Matrix2cd eigvecs = eigs.eigenvectors();
      Vector2cd eigvals = eigs.eigenvalues();
      if ( TD(projector(psiPerp), projector(eigvecs.col(0))) < .1 ) lpsi_psiPerp = real(eigvals(0))*dt;
      else lpsi_psiPerp = real(eigvals(1))*dt;
      if (lpsi_psiPerp >= 0.) psi_psiPerp += lpsi_psiPerp;
      else if (NPerp != 0) psiPerp_psi -= lpsi_psiPerp*(double)Npsi/(double)NPerp;
      eigs.compute(RPerp);
      eigvecs = eigs.eigenvectors();
      eigvals = eigs.eigenvalues();
      if ( TD(projector(psi), projector(eigvecs.col(0))) < .1 ) lpsiPerp_psi = real(eigvals(0))*dt;
      else lpsiPerp_psi = real(eigvals(1))*dt;
      if (lpsiPerp_psi >= 0.) psiPerp_psi += lpsiPerp_psi;
      else if (Npsi != 0) psi_psiPerp -= lpsiPerp_psi*(double)NPerp/(double)Npsi;

      if (Ncopy == 0) {
        out << lpsi_psiPerp << " " << lpsiPerp_psi << endl;
        out << psi_psiPerp << " " << psiPerp_psi << endl;
      }

      for (int j = 0; j < Npsi_old; ++j)
        jump(psi_psiPerp, Npsi, NPerp);
      for (int j = 0; j < NPerp_old; ++j)
        jump(psiPerp_psi, NPerp, Npsi);

      Nstep++;
    }
  }

  // Not printed (and to be changed from the .py): obs and coh for rho and error
  VectorXd m_obs(Nsteps), m_coh(Nsteps), err_obs(Nsteps), err_coh(Nsteps), m_err(Nsteps);
  out.close();
  out.open("driven_obs.txt");
  for (int i = 0; i < Nsteps; ++i) {
    m_obs[i] = observables.col(i).mean();
    m_coh[i] = coherences.col(i).mean();
    err_obs[i] = sqrt((observables.col(i).array() - m_obs[i]).square().sum() / (observables.col(i).size() - 1));
    err_coh[i] = sqrt((coherences.col(i).array() - m_coh[i]).square().sum() / (coherences.col(i).size() - 1));
    m_err[i] = errors.col(i).mean();
    out << m_obs[i] << " " << err_obs[i] << " " << m_coh[i] << " " << err_coh[i] << " " << m_err[i] << endl;
  }

  return 0;
}