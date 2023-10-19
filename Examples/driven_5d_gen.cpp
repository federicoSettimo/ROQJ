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

MatrixXcd projector (const VectorXcd &psi) {return psi*psi.adjoint();}

// No rev
/*double gamma_m (double t) {return 1.;}
double gamma_p (double t) {return .5*exp(-t)*cos(3.*t)+.9;}
double beta (double t) {
  double gp = gamma_p(t), gm = gamma_m(t);
  return .5*min(gp-.5*gm, gm-.5*gp);
}*/
// Rev from pm
/*double gamma_m (double t) {return .5*(cos(t)+1.);}
double gamma_p (double t) {return .5*(-cos(t)+1.);}
double beta (double t) {
  double gp = gamma_p(t), gm = gamma_m(t);
  return .5*min(gp,gm);
}*/
double gamma_m (double t) {return 1.;}
double gamma_p (double t) {return .5*(-cos(t)+1.);}
//double gamma_p (double t) {return exp(-.25*t);}
double beta (double t) {
  double gp = gamma_p(t), gm = gamma_m(t);
  //return t > 1. ? .5*min(gp,gm) : 0.;
  return .5*min(gp,gm);
}
/*double gamma_m (double t) {return 1.;}
double gamma_p (double t) {return 0.;}
double beta (double t) {return .2;}*/
/*double gamma_m (double t) {return 1.;}
double gamma_p (double t) {return .5*(tanh(.5*t)+1);}
double beta (double t) {return .5*gamma_p(t);}*/

Matrix2cd H (double t) {
  return beta(t)*sigma_y;
}

Matrix2cd J (const MatrixXcd &rho, double t) {
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


double tmin = 0., tmax = 10., dt = .001, threshold = 1e-3;
int Nensemble = 10000, Npsi = Nensemble, Nbar = 0, Nbarperp = 0, Np = 0, Nm = 0;
int Npsi_old, Nbar_old, Nbarperp_old, Np_old, Nm_old;
vector<Matrix2cd> exact_dyn((int)(tmax/dt)+1);

void jump (double p1, double p2, int &N, int &N1, int &N2) {
  double z = (double)rand()/(double)RAND_MAX;
  if (z < p1) {N--; N1++;}
  else if (z < p1+p2) {N--; N2++;}
  return;
}

void jump3 (double p1, double p2, double p3, int &N, int &N1, int &N2, int &N3) {
  double z = (double)rand()/(double)RAND_MAX;
  if (z < p1) {N--; N1++;}
  else if (z < p1+p2) {N--; N2++;}
  else if (z < p1+p2+p3) {N--; N3++;}
  return;
}

double rev_jump (double l, int N, int N_target) {
  if (l > 0. || N == 0.) return 0.;
  return -l*(double)N_target/(double)N;
}

Matrix2cd L (const Matrix2cd &rho, double t) {return -I*comm(H(t),rho) + J(rho,t) - .5*anticomm(Gamma(t),rho);}

// Runge-Kutta 4
Matrix2cd RK4 (const Matrix2cd &rho, double t) {
  Matrix2cd k1 = dt * L(rho, t), k2 = dt * L(rho + .5*k1, t + .5*dt), k3 = dt * L(rho + .5*k2, t + .5*dt), k4 = dt * L(rho + k3, t + dt);
  return (1./6.) * (k1 + 2.*k2 + 2.*k3 + k4);
}

// Calculates the error for a particular choice of mu
double get_max_err (const Vector2cd &initialState, double mu);

// Phis
Vector2cd Phi_psi (const Vector2cd &psi, double mu, double t);
Vector2cd Phi_bar (double mu, double t);
Vector2cd Phi_barperp (double mu, double t);
Vector2cd Phi_plus (double mu, double t);
Vector2cd Phi_min (double mu, double t);

int main () {
  srand(time(NULL));
  ofstream out;
  out.open("driven.txt");

  Vector2cd initialState;
  initialState << .7, .8;
  initialState.normalize();

  // Calculating the exact dynamics (with the large timestep)
  exact_dyn[0] = projector(initialState);
  int i = 1;
  for (double t = 0; t < tmax; t += dt) {
    exact_dyn[i] = exact_dyn[i-1] + RK4(exact_dyn[i-1],t);
    ++i;
  }

  // Determining mu that minimizes the error
  cout << "Determining the target states which minimize the average error...\n";
  double best_mu = 0., min_error = 0.;
  double mu_min = -1., mu_max = 1., dmu = 0.01;
  //double mu_min = -.79, mu_max = -.77, dmu = 0.01;
  for (double mu = mu_min; mu < mu_max + dmu; mu += dmu) {
    double err = get_max_err(initialState, mu);
    if (err < min_error || mu == mu_min) {
      min_error = err;
      best_mu = mu;
    }
    cout << "\t" << mu << " " << err << endl;
  }
  double mu = best_mu;
  //dt = dt/50.; // Using a smaller timestep for the actual simulations
  cout << "Best state: " << sqrt(1.-mu*mu) << " |1> + " << mu << " |0>\n";
  out << (int)(tmax/dt) << endl << tmax << endl << dt << endl;
  out << mu << endl;

  Vector2cd psi = initialState.normalized(), psibar, psibarperp, psiplus, psimin;
  Matrix2cd rho = projector(psi), exact = projector(psi);
  psibar = sqrt(1.-mu*mu)*excited_state + mu*ground_state;
  psibarperp = sqrt(1.-mu*mu)*ground_state - mu*excited_state;
  psiplus = (psibar + psibarperp)/sqrt(2);
  psimin = (psibar - psibarperp)/sqrt(2);

  Npsi = Nensemble; Nbar = 0; Nbarperp = 0; Np = 0; Nm = 0;

  for (double t = 0.; t < tmax; t += dt) {
    Npsi_old = Npsi; Nbar_old = Nbar; Nbarperp_old = Nbarperp; Np_old = Np; Nm_old = Nm;
    double fpsi = (double)Npsi/(double)Nensemble, fp = (double)Np/(double)Nensemble, fm = (double)Nm/(double)Nensemble, fbar = (double)Nbar/(double)Nensemble, fbarperp = (double)Nbarperp/(double)Nensemble;
    rho = fpsi*projector(psi) + fbar*projector(psibar) + fbarperp*projector(psibarperp) + fp*projector(psiplus) + fm*projector(psimin);
    double gp = gamma_p(t), gm = gamma_m(t), b = beta(t);

    out << observable(rho) << " " << observable(exact) << " " << observable(projector(psi)) << endl;
    out << coherence(rho) << " " << coherence(exact) << " " << coherence(projector(psi)) << endl;
    out << TD(rho, exact) << endl;
    out << fpsi << " " << fbar << " " << fbarperp << " " << fp << " " << fm << endl;
    out << gp << " " << gm << " " << b << endl;
    cout << TD(rho, exact) << endl;

    // Updating exact
    exact += RK4(exact,t);

    // Updating the occupation probabilities: first getting R
    Vector2cd phi = Phi_psi(psi, mu, t), phi_bar = Phi_bar(mu, t), phi_barperp = Phi_barperp(mu, t), phi_plus = Phi_plus(mu, t), phi_min = Phi_min(mu, t);
    Matrix2cd Rpsi = J(projector(psi),t) + .5*(psi*phi.adjoint() + phi*psi.adjoint());
    Matrix2cd Rbar = J(projector(psibar),t) + .5*(psibar*phi_bar.adjoint() + phi_bar*psibar.adjoint());
    Matrix2cd Rbarperp = J(projector(psibarperp),t) + .5*(psibarperp*phi_barperp.adjoint() + phi_barperp*psibarperp.adjoint());
    Matrix2cd Rplus = J(projector(psiplus),t) + .5*(psiplus*phi_plus.adjoint() + phi_plus*psiplus.adjoint());
    Matrix2cd Rmin = J(projector(psimin),t) + .5*(psimin*phi_min.adjoint() + phi_min*psimin.adjoint());
    // Then, calculating the rates
    double lpsi_bar, lpsi_barperp; // From psi
    ComplexEigenSolver<Matrix2cd> eigs; eigs.compute(Rpsi); Vector2cd eigval = eigs.eigenvalues(); Matrix2cd eigvec = eigs.eigenvectors();
    if ( (eigvec.col(0) - psibar).norm() < .1 ) {lpsi_bar = real(eigval(0))*dt; lpsi_barperp = real(eigval(1))*dt;}
    else {lpsi_bar = real(eigval(1))*dt; lpsi_barperp = real(eigval(0))*dt;}
    double lbar_plus, lbar_min; // From psibar
    eigs.compute(Rbar); eigval = eigs.eigenvalues(); eigvec = eigs.eigenvectors();
    if ( (eigvec.col(0) - psiplus).norm() < .1 ) {lbar_plus = real(eigval(0))*dt; lbar_min = real(eigval(1))*dt;}
    else {lbar_plus = real(eigval(1))*dt; lbar_min = real(eigval(0))*dt;}
    double lbarperp_plus, lbarperp_min; // From psibarperp
    eigs.compute(Rbarperp); eigval = eigs.eigenvalues(); eigvec = eigs.eigenvectors();
    if ( (eigvec.col(0) - psiplus).norm() < .1 ) {lbarperp_plus = real(eigval(0))*dt; lbarperp_min = real(eigval(1))*dt;}
    else {lbarperp_plus = real(eigval(1))*dt; lbarperp_min = real(eigval(0))*dt;}
    double lplus_bar, lplus_barperp; // From psiplus
    eigs.compute(Rplus); eigval = eigs.eigenvalues(); eigvec = eigs.eigenvectors();
    if ( (eigvec.col(0) - psibar).norm() < .1 ) {lplus_bar = real(eigval(0))*dt; lplus_barperp = real(eigval(1))*dt;}
    else {lplus_bar = real(eigval(1))*dt; lplus_barperp = real(eigval(0))*dt;}
    double lmin_bar, lmin_barperp; // From psimin
    eigs.compute(Rmin); eigval = eigs.eigenvalues(); eigvec = eigs.eigenvectors();
    if ( (eigvec.col(0) - psibar).norm() < .1 ) {lmin_bar = real(eigval(0))*dt; lmin_barperp = real(eigval(1))*dt;}
    else {lmin_bar = real(eigval(1))*dt; lmin_barperp = real(eigval(0))*dt;}

    // Positive or negative rates?
    double psi_bar = 0., psi_barperp = 0., bar_plus = 0., bar_min = 0., barperp_plus = 0., barperp_min = 0., plus_bar = 0., plus_barperp = 0., min_bar = 0., min_barperp = 0., bar_psi = 0., barperp_psi = 0.;
    if (lpsi_bar >= 0.) psi_bar += lpsi_bar; // psi
    else {bar_psi += rev_jump(lpsi_bar, Nbar, Npsi); barperp_psi += rev_jump(lpsi_bar, Nbarperp, Npsi);}
    if (lpsi_barperp >= 0.) psi_barperp += lpsi_barperp;
    else {barperp_psi += rev_jump(lpsi_barperp, Nbar, Npsi); barperp_psi += rev_jump(lpsi_barperp, Nbarperp, Npsi);}
    if (lbar_plus >= 0.) bar_plus += lbar_plus; // psibar
    else {psi_bar += rev_jump(lbar_plus, Npsi, Nbar); plus_bar += rev_jump(lbar_plus, Np, Nbar); min_bar += rev_jump(lbar_plus, Nm, Nbar);}
    if (lbar_min >= 0.) bar_min += lbar_min;
    else {psi_bar += rev_jump(lbar_min, Npsi, Nbar); plus_bar += rev_jump(lbar_min, Np, Nbar); min_bar += rev_jump(lbar_min, Nm, Nbar);}
    if (lbarperp_plus >= 0.) barperp_plus += lbarperp_plus; // psibarperp
    else {psi_barperp += rev_jump(lbarperp_plus, Npsi, Nbarperp); plus_barperp += rev_jump(lbarperp_plus, Np, Nbarperp); min_barperp += rev_jump(lbarperp_plus, Nm, Nbarperp);}
    if (lbarperp_min >= 0.) barperp_min += lbarperp_min;
    else {psi_barperp += rev_jump(lbarperp_min, Npsi, Nbarperp); plus_barperp += rev_jump(lbarperp_min, Np, Nbarperp); min_barperp += rev_jump(lbarperp_min, Nm, Nbarperp);}
    if (lplus_bar >= 0.) plus_bar += lplus_bar; // psiplus
    else {bar_plus += rev_jump(lplus_bar, Nbar, Np); barperp_plus += rev_jump(lplus_bar, Nbarperp, Np);}
    if (lplus_barperp >= 0.) plus_barperp += lplus_barperp;
    else {bar_plus += rev_jump(lplus_barperp, Nbar, Np); barperp_plus += rev_jump(lplus_barperp, Nbarperp, Np);}
    if (lmin_bar >= 0.) min_bar += lmin_bar; // psimin
    else {bar_min += rev_jump(lmin_bar, Nbar, Nm); barperp_min += rev_jump(lmin_bar, Nbarperp, Nm);}
    if (lmin_barperp >= 0.) min_barperp += lmin_barperp;
    else {bar_min += rev_jump(lmin_barperp, Nbar, Nm); barperp_min += rev_jump(lmin_barperp, Nbarperp, Nm);}

    // Now doing the jumps
    for (int i = 0; i < Npsi_old; ++i)
      jump(psi_bar, psi_barperp, Npsi, Nbar, Nbarperp);
    for (int i = 0; i < Nbar_old; ++i) 
      jump3(bar_plus, bar_min, bar_psi, Nbar, Np, Nm, Npsi);
    for (int i = 0; i < Nbarperp_old; ++i) 
      jump3(barperp_plus, barperp_min, barperp_psi, Nbarperp, Np, Nm, Npsi);
    for (int i = 0; i < Np_old; ++i)
      jump (plus_bar, plus_barperp, Np, Nbar, Nbarperp);
    for (int i = 0; i < Nm_old; ++i)
      jump (min_bar, min_barperp, Nm, Nbar, Nbarperp);

    // Deterministically evolve psi
    Matrix2cd K = H(t) - .5*I*Gamma(t);
    psi -= dt*(I*K*psi + .5*phi);
    psi.normalize();

    out << lpsi_bar << " " << lpsi_barperp << " " << lbarperp_plus << " " << lbarperp_min << " " << lbar_plus << " " << bar_min << " " << lplus_bar << " " << lplus_barperp << " " << lmin_bar << " " << lmin_barperp << endl;
    out << psi_bar << " " << psi_barperp << " " << barperp_plus << " " << barperp_min << " " << bar_plus << " " << bar_min << " " << plus_bar << " " << plus_barperp << " " << min_bar << " " << min_barperp << " " << bar_psi << " " << barperp_psi << endl;
  }

  return 0;
}


double get_max_err (const Vector2cd &initialState, double mu) {
  Vector2cd psi = initialState.normalized(), psibar, psibarperp, psiplus, psimin;
  Matrix2cd rho = projector(psi), exact = projector(psi);
  psibar << sqrt(1.-mu*mu), mu; psibar.normalize();
  psibarperp << -mu, sqrt(1.-mu*mu); psibarperp.normalize();
  psiplus = (psibar + psibarperp)/sqrt(2);
  psimin = (psibar - psibarperp)/sqrt(2);

  Npsi = Nensemble; Nbar = 0; Nbarperp = 0; Np = 0; Nm = 0;

  double max_err = 0.;
  int i = 0;
  for (double t = 0.; t < tmax; t += dt) {
    Npsi_old = Npsi; Nbar_old = Nbar; Nbarperp_old = Nbarperp; Np_old = Np; Nm_old = Nm;
    double fpsi = (double)Npsi/(double)Nensemble, fp = (double)Np/(double)Nensemble, fm = (double)Nm/(double)Nensemble, fbar = (double)Nbar/(double)Nensemble, fbarperp = (double)Nbarperp/(double)Nensemble;
    rho = fpsi*projector(psi) + fbar*projector(psibar) + fbarperp*projector(psibarperp) + fp*projector(psiplus) + fm*projector(psimin);

    double err = TD(rho, exact_dyn[i]);
    if (isnan(err) || isnan(observable(rho)) || isnan(coherence(rho))) return 1.;
    if (err > max_err)
      max_err = err;
    ++i;

    // Updating the occupation probabilities: first getting R
    Vector2cd phi = Phi_psi(psi, mu, t), phi_bar = Phi_bar(mu, t), phi_barperp = Phi_barperp(mu, t), phi_plus = Phi_plus(mu, t), phi_min = Phi_min(mu, t);
    Matrix2cd Rpsi = J(projector(psi),t) + .5*(psi*phi.adjoint() + phi*psi.adjoint());
    Matrix2cd Rbar = J(projector(psibar),t) + .5*(psibar*phi_bar.adjoint() + phi_bar*psibar.adjoint());
    Matrix2cd Rbarperp = J(projector(psibarperp),t) + .5*(psibarperp*phi_barperp.adjoint() + phi_barperp*psibarperp.adjoint());
    Matrix2cd Rplus = J(projector(psiplus),t) + .5*(psiplus*phi_plus.adjoint() + phi_plus*psiplus.adjoint());
    Matrix2cd Rmin = J(projector(psimin),t) + .5*(psimin*phi_min.adjoint() + phi_min*psimin.adjoint());
    // Then, calculating the rates
    double lpsi_bar, lpsi_barperp; // From psi
    ComplexEigenSolver<Matrix2cd> eigs; eigs.compute(Rpsi); Vector2cd eigval = eigs.eigenvalues(); Matrix2cd eigvec = eigs.eigenvectors();
    if ( (eigvec.col(0) - psibar).norm() < .1 ) {lpsi_bar = real(eigval(0))*dt; lpsi_barperp = real(eigval(1))*dt;}
    else {lpsi_bar = real(eigval(1))*dt; lpsi_barperp = real(eigval(0))*dt;}
    double lbar_plus, lbar_min; // From psibar
    eigs.compute(Rbar); eigval = eigs.eigenvalues(); eigvec = eigs.eigenvectors();
    if ( (eigvec.col(0) - psiplus).norm() < .1 ) {lbar_plus = real(eigval(0))*dt; lbar_min = real(eigval(1))*dt;}
    else {lbar_plus = real(eigval(1))*dt; lbar_min = real(eigval(0))*dt;}
    double lbarperp_plus, lbarperp_min; // From psibarperp
    eigs.compute(Rbarperp); eigval = eigs.eigenvalues(); eigvec = eigs.eigenvectors();
    if ( (eigvec.col(0) - psiplus).norm() < .1 ) {lbarperp_plus = real(eigval(0))*dt; lbarperp_min = real(eigval(1))*dt;}
    else {lbarperp_plus = real(eigval(1))*dt; lbarperp_min = real(eigval(0))*dt;}
    double lplus_bar, lplus_barperp; // From psiplus
    eigs.compute(Rplus); eigval = eigs.eigenvalues(); eigvec = eigs.eigenvectors();
    if ( (eigvec.col(0) - psibar).norm() < .1 ) {lplus_bar = real(eigval(0))*dt; lplus_barperp = real(eigval(1))*dt;}
    else {lplus_bar = real(eigval(1))*dt; lplus_barperp = real(eigval(0))*dt;}
    double lmin_bar, lmin_barperp; // From psimin
    eigs.compute(Rmin); eigval = eigs.eigenvalues(); eigvec = eigs.eigenvectors();
    if ( (eigvec.col(0) - psibar).norm() < .1 ) {lmin_bar = real(eigval(0))*dt; lmin_barperp = real(eigval(1))*dt;}
    else {lmin_bar = real(eigval(1))*dt; lmin_barperp = real(eigval(0))*dt;}

    // Positive or negative rates?
    double psi_bar = 0., psi_barperp = 0., bar_plus = 0., bar_min = 0., barperp_plus = 0., barperp_min = 0., plus_bar = 0., plus_barperp = 0., min_bar = 0., min_barperp = 0., bar_psi = 0., barperp_psi = 0.;
    if (lpsi_bar >= 0.) psi_bar += lpsi_bar; // psi
    else {bar_psi += rev_jump(lpsi_bar, Nbar, Npsi); barperp_psi += rev_jump(lpsi_bar, Nbarperp, Npsi);}
    if (lpsi_barperp >= 0.) psi_barperp += lpsi_barperp;
    else {barperp_psi += rev_jump(lpsi_barperp, Nbar, Npsi); barperp_psi += rev_jump(lpsi_barperp, Nbarperp, Npsi);}
    if (lbar_plus >= 0.) bar_plus += lbar_plus; // psibar
    else {psi_bar += rev_jump(lbar_plus, Npsi, Nbar); plus_bar += rev_jump(lbar_plus, Np, Nbar); min_bar += rev_jump(lbar_plus, Nm, Nbar);}
    if (lbar_min >= 0.) bar_min += lbar_min;
    else {psi_bar += rev_jump(lbar_min, Npsi, Nbar); plus_bar += rev_jump(lbar_min, Np, Nbar); min_bar += rev_jump(lbar_min, Nm, Nbar);}
    if (lbarperp_plus >= 0.) barperp_plus += lbarperp_plus; // psibarperp
    else {psi_barperp += rev_jump(lbarperp_plus, Npsi, Nbarperp); plus_barperp += rev_jump(lbarperp_plus, Np, Nbarperp); min_barperp += rev_jump(lbarperp_plus, Nm, Nbarperp);}
    if (lbarperp_min >= 0.) barperp_min += lbarperp_min;
    else {psi_barperp += rev_jump(lbarperp_min, Npsi, Nbarperp); plus_barperp += rev_jump(lbarperp_min, Np, Nbarperp); min_barperp += rev_jump(lbarperp_min, Nm, Nbarperp);}
    if (lplus_bar >= 0.) plus_bar += lplus_bar; // psiplus
    else {bar_plus += rev_jump(lplus_bar, Nbar, Np); barperp_plus += rev_jump(lplus_bar, Nbarperp, Np);}
    if (lplus_barperp >= 0.) plus_barperp += lplus_barperp;
    else {bar_plus += rev_jump(lplus_barperp, Nbar, Np); barperp_plus += rev_jump(lplus_barperp, Nbarperp, Np);}
    if (lmin_bar >= 0.) min_bar += lmin_bar; // psimin
    else {bar_min += rev_jump(lmin_bar, Nbar, Nm); barperp_min += rev_jump(lmin_bar, Nbarperp, Nm);}
    if (lmin_barperp >= 0.) min_barperp += lmin_barperp;
    else {bar_min += rev_jump(lmin_barperp, Nbar, Nm); barperp_min += rev_jump(lmin_barperp, Nbarperp, Nm);}

    // Now doing the jumps
    for (int i = 0; i < Npsi_old; ++i)
      jump(psi_bar, psi_barperp, Npsi, Nbar, Nbarperp);
    for (int i = 0; i < Nbar_old; ++i) 
      jump3(bar_plus, bar_min, bar_psi, Nbar, Np, Nm, Npsi);
    for (int i = 0; i < Nbarperp_old; ++i) 
      jump3(barperp_plus, barperp_min, barperp_psi, Nbarperp, Np, Nm, Npsi);
    for (int i = 0; i < Np_old; ++i)
      jump (plus_bar, plus_barperp, Np, Nbar, Nbarperp);
    for (int i = 0; i < Nm_old; ++i)
      jump (min_bar, min_barperp, Nm, Nbar, Nbarperp);

    // Deterministically evolve psi
    Matrix2cd K = H(t) - .5*I*Gamma(t);
    psi -= dt*(I*K*psi + .5*phi);
    psi.normalize();
  }

  return max_err;
}

Vector2cd Phi_psi (const Vector2cd &psi, double mu, double t) {
  double a = real(psi(1)), gp = gamma_p(t), gm = gamma_m(t), b = beta(t);
  //double ph1 = sqrt(1.-a*a);
  double ph1 = (sqrt(1 - pow(a,2))*gm*pow(mu,2) - 3*sqrt(1 - pow(a,2))*gm*pow(mu,4) + 4*sqrt(1 - pow(a,2))*gm*pow(mu,6) - 2*sqrt(1 - pow(a,2))*gm*pow(mu,8) - 
     2*pow(a,9)*(gm + gp)*mu*sqrt(1 - pow(mu,2))*(-1 + 2*pow(mu,2)) + 4*pow(a,7)*(2*gm + gp)*mu*sqrt(1 - pow(mu,2))*(-1 + 2*pow(mu,2)) + 
     2*a*gm*mu*sqrt(1 - pow(mu,2))*(1 - 3*pow(mu,2) + 3*pow(mu,4) - 2*pow(mu,6)) + 
     4*pow(a,3)*gm*mu*sqrt(1 - pow(mu,2))*(-2 + 5*pow(mu,2) - 3*pow(mu,4) + 2*pow(mu,6)) + 
     pow(a,8)*sqrt(1 - pow(a,2))*(gp + 4*gp*pow(mu,2) - 4*gp*pow(mu,4) - gm*pow(1 - 2*pow(mu,2),2)) - 
     2*pow(a,5)*mu*sqrt(1 - pow(mu,2))*(-1 + 2*pow(mu,2))*(gp*(1 - pow(mu,2) + pow(mu,4)) + gm*(6 - pow(mu,2) + pow(mu,4))) + 
     pow(a,6)*sqrt(1 - pow(a,2))*(gp*(-1 - 6*pow(mu,2) + 6*pow(mu,4)) + gm*(3 - 14*pow(mu,2) + 14*pow(mu,4))) + 
     pow(a,4)*sqrt(1 - pow(a,2))*(gp*pow(mu,2)*(5 - 9*pow(mu,2) + 8*pow(mu,4) - 4*pow(mu,6)) + 
        gm*(-3 + 17*pow(mu,2) - 21*pow(mu,4) + 8*pow(mu,6) - 4*pow(mu,8))) + 
     pow(a,2)*sqrt(1 - pow(a,2))*(-2*gp*pow(mu,4)*pow(-1 + pow(mu,2),2) + gm*(1 - 8*pow(mu,2) + 14*pow(mu,4) - 12*pow(mu,6) + 6*pow(mu,8))) - 
     sqrt(-(pow(pow(-1 + pow(a,2),2)*gm + pow(a,4)*gp,2)*(4*pow(a,9)*sqrt(1 - pow(a,2))*mu*(1 - 2*pow(mu,2))*sqrt(1 - pow(mu,2)) - 
           pow(mu,4)*pow(1 - 3*pow(mu,2) + 2*pow(mu,4),2) - 
           4*pow(a,5)*sqrt(1 - pow(a,2))*mu*sqrt(1 - pow(mu,2))*pow(-1 + 2*pow(mu,2),3)*(3 - 4*pow(mu,2) + 4*pow(mu,4)) + 
           pow(a,10)*(1 - 8*pow(mu,2) + 8*pow(mu,4)) + 4*pow(a,7)*sqrt(1 - pow(a,2))*mu*sqrt(1 - pow(mu,2))*
            (-3 + 14*pow(mu,2) - 24*pow(mu,4) + 16*pow(mu,6)) + pow(a,8)*(-3 + 36*pow(mu,2) - 100*pow(mu,4) + 128*pow(mu,6) - 64*pow(mu,8)) - 
           4*a*sqrt(1 - pow(a,2))*pow(mu,3)*sqrt(1 - pow(mu,2))*(1 - 6*pow(mu,2) + 14*pow(mu,4) - 15*pow(mu,6) + 6*pow(mu,8)) + 
           4*pow(a,3)*sqrt(1 - pow(a,2))*mu*sqrt(1 - pow(mu,2))*(-1 + 11*pow(mu,2) - 43*pow(mu,4) + 82*pow(mu,6) - 80*pow(mu,8) + 32*pow(mu,10)) + 
           pow(a,2)*pow(mu,2)*(-6 + 47*pow(mu,2) - 154*pow(mu,4) + 257*pow(mu,6) - 216*pow(mu,8) + 72*pow(mu,10)) + 
           pow(a,4)*(-1 + 32*pow(mu,2) - 184*pow(mu,4) + 496*pow(mu,6) - 728*pow(mu,8) + 576*pow(mu,10) - 192*pow(mu,12)) + 
           pow(a,6)*(3 - 54*pow(mu,2) + 230*pow(mu,4) - 480*pow(mu,6) + 560*pow(mu,8) - 384*pow(mu,10) + 128*pow(mu,12))))))/
   (2.*pow(pow(a,2) - pow(a,4) - pow(mu,2) + pow(mu,4),2));
  return ph1*excited_state +
    ((-2*(-1 + pow(a,2))*gm*mu*sqrt(1 - pow(mu,2)) - 2*pow(a,2)*gp*mu*sqrt(1 - pow(mu,2)) - 2*sqrt(1 - pow(a,2))*mu*sqrt(1 - pow(mu,2))*ph1 + 
     a*(ph1 - 2*pow(mu,2)*ph1))/(-sqrt(1 - pow(a,2)) + 2*sqrt(1 - pow(a,2))*pow(mu,2) - 2*a*mu*sqrt(1 - pow(mu,2))))*ground_state;
}
Vector2cd Phi_bar (double mu, double t) {
  double gp = gamma_p(t), gm = gamma_m(t), b = beta(t);
  return (-2*b*mu + sqrt(1 - pow(mu,2))*(gm - 4*gm*pow(mu,2) + 2*gm*pow(mu,4) + 2*gp*pow(mu,4)))*excited_state +
    (2*b*sqrt(1 - pow(mu,2)) + 2*gm*mu*pow(-1 + pow(mu,2),2) + gp*mu*(-1 + 2*pow(mu,4)))*ground_state;
}
Vector2cd Phi_barperp (double mu, double t) {
  double gp = gamma_p(t), gm = gamma_m(t), b = beta(t);
  return (gm*(mu - 2*pow(mu,5)) - 2*(b*sqrt(1 - pow(mu,2)) + gp*mu*pow(-1 + pow(mu,2),2)))*excited_state +
    (-2*b*mu + sqrt(1 - pow(mu,2))*(gp - 4*gp*pow(mu,2) + 2*gm*pow(mu,4) + 2*gp*pow(mu,4)))*ground_state;
}
Vector2cd Phi_plus (double mu, double t) {
  double gp = gamma_p(t), gm = gamma_m(t), b = beta(t);
  return ((-4*b*(1 + 2*mu*sqrt(1 - pow(mu,2))) + (-1 + 2*pow(mu,2))*(-(gp*(1 + 4*pow(mu,2) - 4*pow(mu,4) + 4*mu*sqrt(1 - pow(mu,2)))) + 
        gm*(1 - 4*pow(mu,2) + 4*pow(mu,4) + 4*mu*sqrt(1 - pow(mu,2)))))/(2.*sqrt(2)*(mu + sqrt(1 - pow(mu,2)))))*excited_state +
    ((4*b*(-mu + sqrt(1 - pow(mu,2))) + gm*(-3*mu + 8*pow(mu,3) - 4*pow(mu,5) + sqrt(1 - pow(mu,2)) - 4*pow(mu,4)*sqrt(1 - pow(mu,2))) - 
     gp*(-3*mu + 4*pow(mu,5) + sqrt(1 - pow(mu,2)) - 8*pow(mu,2)*sqrt(1 - pow(mu,2)) + 4*pow(mu,4)*sqrt(1 - pow(mu,2))))/(2.*sqrt(2)))*ground_state;
}
Vector2cd Phi_min (double mu, double t) {
  double gp = gamma_p(t), gm = gamma_m(t), b = beta(t);
  return ((b*(4 - 8*mu*sqrt(1 - pow(mu,2))) + (-1 + 2*pow(mu,2))*(gm*(1 - 4*pow(mu,2) + 4*pow(mu,4) - 4*mu*sqrt(1 - pow(mu,2))) + 
        gp*(-1 - 4*pow(mu,2) + 4*pow(mu,4) + 4*mu*sqrt(1 - pow(mu,2)))))/(2.*sqrt(2)*(-mu + sqrt(1 - pow(mu,2)))))*excited_state +
    ((4*b*(mu + sqrt(1 - pow(mu,2))) + gm*(-1 + 2*pow(mu,2))*(3*mu - 2*pow(mu,3) + sqrt(1 - pow(mu,2)) + 2*pow(mu,2)*sqrt(1 - pow(mu,2))) + 
     gp*(3*mu - 4*pow(mu,5) + sqrt(1 - pow(mu,2)) - 8*pow(mu,2)*sqrt(1 - pow(mu,2)) + 4*pow(mu,4)*sqrt(1 - pow(mu,2))))/(2.*sqrt(2)))*ground_state;
}