#include "../roqj_state.h"

using namespace std;
using namespace Eigen;

double tmin = 0., tmax = 10., dt = 0.01, threshold = 1e-3;
vector<Vector2cd> phi((int)((tmax-tmin)/dt+1)), psi((int)((tmax-tmin)/dt+1)), phi_perp((int)((tmax-tmin)/dt+1));

double sgn (double x) {return x >= 0. ? 1. : -1.;}

double gamma_m (double t) {return 1;}
double gamma_p (double t) {return 0.;}
double beta (double t) {return 1.;}

MatrixXcd H (double t) {
  return beta(t)*sigma_y;
}

MatrixXcd J (const MatrixXcd &rho, double t) {
  return gamma_m(t)*sigma_m*rho*sigma_p + gamma_p(t)*sigma_p*rho*sigma_m;
}

MatrixXcd Gamma (double t) {
  return gamma_m(t)*sigma_p*sigma_m + gamma_p(t)*sigma_m*sigma_p;
}

VectorXcd Phi (const VectorXcd &psi, double t, bool jumped) {return 0.*ground_state;}

Vector2cd Phi_psi (const VectorXcd &psi, double t, double mu) {
  if (t == 0.) return 0.*ground_state;
  double gm = gamma_m(t), b = beta(t), a = real(psi(1));

  // Case mu = 0
  if (abs(mu) <= .00001)
    return (pow(1 - pow(a,2),1.5)*gm)*excited_state + (a*(-1 + pow(a,2))*gm)*ground_state;

  // Case mu = 1
  if (abs(mu) >= .99999)
    return (-((sqrt(-(pow(-1 + pow(a,2),5)*pow(dt,2)*pow(gm,2))) + pow(a,4)*(-2 + dt)*gm*sqrt(-pow(a,2) + mu) - 2*pow(a,2)*(-1 + dt)*gm*sqrt(-pow(a,2) + mu) + 
       dt*gm*sqrt(-pow(a,2) + mu))/(2*pow(a,4)*(-1 + dt) - 2*pow(a,2)*dt)))*excited_state + 
       (-0.5*(pow(a,4)*(-2 + dt)*gm - 2*pow(a,2)*(-1 + dt)*gm + dt*gm + sqrt(-(pow(-1 + pow(a,2),5)*pow(dt,2)*pow(gm,2)))/sqrt(-pow(a,2) + mu))/
    (a*(-(pow(a,2)*(-1 + dt)) + dt)))*ground_state;

   return ((pow(a,9)*sqrt(1 - pow(a,2))*gm*(1 - 2*pow(mu,2)) + 2*gm*pow(mu,3)*sqrt(1 - pow(mu,2)) - 6*gm*pow(mu,5)*sqrt(1 - pow(mu,2)) + 
     8*gm*pow(mu,7)*sqrt(1 - pow(mu,2)) - 4*gm*pow(mu,9)*sqrt(1 - pow(mu,2)) - 12*pow(a,6)*gm*pow(mu,3)*pow(1 - pow(mu,2),1.5) + 
     4*pow(a,8)*gm*pow(mu,3)*pow(1 - pow(mu,2),1.5) + pow(a,7)*sqrt(1 - pow(a,2))*gm*(-3 + 4*pow(mu,2) + 6*pow(mu,4) - 4*pow(mu,6)) + 
     8*pow(a,2)*gm*pow(mu,3)*sqrt(1 - pow(mu,2))*(-1 + 2*pow(mu,2) - 2*pow(mu,4) + pow(mu,6)) - 
     2*pow(a,4)*gm*pow(mu,3)*sqrt(1 - pow(mu,2))*(-7 + 9*pow(mu,2) - 4*pow(mu,4) + 2*pow(mu,6)) + 
     pow(a,5)*sqrt(1 - pow(a,2))*gm*(3 + pow(mu,2) - 21*pow(mu,4) + 14*pow(mu,6)) - 
     pow(a,3)*sqrt(1 - pow(a,2))*gm*(1 + 6*pow(mu,2) - 26*pow(mu,4) + 24*pow(mu,6) - 10*pow(mu,8) + 4*pow(mu,10)) - 
     sqrt(1 - pow(a,2))*sqrt(pow(-1 + pow(a,2),4)*pow(gm,2)*
        (4*pow(a,9)*sqrt(1 - pow(a,2))*mu*(1 - 2*pow(mu,2))*sqrt(1 - pow(mu,2)) - 4*pow(mu,6)*pow(-1 + pow(mu,2),3) - 
          12*a*sqrt(1 - pow(a,2))*pow(mu,5)*pow(1 - pow(mu,2),2.5)*(-1 + 2*pow(mu,2)) + pow(a,10)*(1 - 8*pow(mu,2) + 8*pow(mu,4)) + 
          9*pow(a,2)*pow(mu,4)*pow(-1 + pow(mu,2),2)*(1 - 8*pow(mu,2) + 8*pow(mu,4)) - 
          4*pow(a,7)*sqrt(1 - pow(a,2))*mu*sqrt(1 - pow(mu,2))*(1 + 6*pow(mu,2) - 24*pow(mu,4) + 16*pow(mu,6)) - 
          32*pow(a,5)*sqrt(1 - pow(a,2))*pow(mu,3)*sqrt(1 - pow(mu,2))*(-1 + pow(mu,2) + 6*pow(mu,4) - 10*pow(mu,6) + 4*pow(mu,8)) + 
          4*pow(a,3)*sqrt(1 - pow(a,2))*pow(mu,3)*sqrt(1 - pow(mu,2))*(-1 - 13*pow(mu,2) + 62*pow(mu,4) - 80*pow(mu,6) + 32*pow(mu,8)) + 
          2*pow(a,8)*(-1 + 2*pow(mu,2) + 30*pow(mu,4) - 64*pow(mu,6) + 32*pow(mu,8)) - 
          2*pow(a,4)*pow(mu,2)*(3 - 7*pow(mu,2) - 88*pow(mu,4) + 284*pow(mu,6) - 288*pow(mu,8) + 96*pow(mu,10)) + 
          pow(a,6)*(1 + 10*pow(mu,2) - 90*pow(mu,4) + 32*pow(mu,6) + 304*pow(mu,8) - 384*pow(mu,10) + 128*pow(mu,12)))) + 
     2*sqrt(1 - pow(a,2))*pow(mu,2)*sqrt(pow(-1 + pow(a,2),4)*pow(gm,2)*
        (4*pow(a,9)*sqrt(1 - pow(a,2))*mu*(1 - 2*pow(mu,2))*sqrt(1 - pow(mu,2)) - 4*pow(mu,6)*pow(-1 + pow(mu,2),3) - 
          12*a*sqrt(1 - pow(a,2))*pow(mu,5)*pow(1 - pow(mu,2),2.5)*(-1 + 2*pow(mu,2)) + pow(a,10)*(1 - 8*pow(mu,2) + 8*pow(mu,4)) + 
          9*pow(a,2)*pow(mu,4)*pow(-1 + pow(mu,2),2)*(1 - 8*pow(mu,2) + 8*pow(mu,4)) - 
          4*pow(a,7)*sqrt(1 - pow(a,2))*mu*sqrt(1 - pow(mu,2))*(1 + 6*pow(mu,2) - 24*pow(mu,4) + 16*pow(mu,6)) - 
          32*pow(a,5)*sqrt(1 - pow(a,2))*pow(mu,3)*sqrt(1 - pow(mu,2))*(-1 + pow(mu,2) + 6*pow(mu,4) - 10*pow(mu,6) + 4*pow(mu,8)) + 
          4*pow(a,3)*sqrt(1 - pow(a,2))*pow(mu,3)*sqrt(1 - pow(mu,2))*(-1 - 13*pow(mu,2) + 62*pow(mu,4) - 80*pow(mu,6) + 32*pow(mu,8)) + 
          2*pow(a,8)*(-1 + 2*pow(mu,2) + 30*pow(mu,4) - 64*pow(mu,6) + 32*pow(mu,8)) - 
          2*pow(a,4)*pow(mu,2)*(3 - 7*pow(mu,2) - 88*pow(mu,4) + 284*pow(mu,6) - 288*pow(mu,8) + 96*pow(mu,10)) + 
          pow(a,6)*(1 + 10*pow(mu,2) - 90*pow(mu,4) + 32*pow(mu,6) + 304*pow(mu,8) - 384*pow(mu,10) + 128*pow(mu,12)))) + 
     a*mu*(sqrt(1 - pow(a,2))*gm*mu*(3 - 11*pow(mu,2) + 14*pow(mu,4) - 10*pow(mu,6) + 4*pow(mu,8)) - 
        2*sqrt(1 - pow(mu,2))*sqrt(pow(-1 + pow(a,2),4)*pow(gm,2)*
           (4*pow(a,9)*sqrt(1 - pow(a,2))*mu*(1 - 2*pow(mu,2))*sqrt(1 - pow(mu,2)) - 4*pow(mu,6)*pow(-1 + pow(mu,2),3) - 
             12*a*sqrt(1 - pow(a,2))*pow(mu,5)*pow(1 - pow(mu,2),2.5)*(-1 + 2*pow(mu,2)) + pow(a,10)*(1 - 8*pow(mu,2) + 8*pow(mu,4)) + 
             9*pow(a,2)*pow(mu,4)*pow(-1 + pow(mu,2),2)*(1 - 8*pow(mu,2) + 8*pow(mu,4)) - 
             4*pow(a,7)*sqrt(1 - pow(a,2))*mu*sqrt(1 - pow(mu,2))*(1 + 6*pow(mu,2) - 24*pow(mu,4) + 16*pow(mu,6)) - 
             32*pow(a,5)*sqrt(1 - pow(a,2))*pow(mu,3)*sqrt(1 - pow(mu,2))*(-1 + pow(mu,2) + 6*pow(mu,4) - 10*pow(mu,6) + 4*pow(mu,8)) + 
             4*pow(a,3)*sqrt(1 - pow(a,2))*pow(mu,3)*sqrt(1 - pow(mu,2))*(-1 - 13*pow(mu,2) + 62*pow(mu,4) - 80*pow(mu,6) + 32*pow(mu,8)) + 
             2*pow(a,8)*(-1 + 2*pow(mu,2) + 30*pow(mu,4) - 64*pow(mu,6) + 32*pow(mu,8)) - 
             2*pow(a,4)*pow(mu,2)*(3 - 7*pow(mu,2) - 88*pow(mu,4) + 284*pow(mu,6) - 288*pow(mu,8) + 96*pow(mu,10)) + 
             pow(a,6)*(1 + 10*pow(mu,2) - 90*pow(mu,4) + 32*pow(mu,6) + 304*pow(mu,8) - 384*pow(mu,10) + 128*pow(mu,12))))))/
   (2.*pow(pow(a,2) - pow(a,4) - pow(mu,2) + pow(mu,4),2)*(2*sqrt(1 - pow(a,2))*mu*sqrt(1 - pow(mu,2)) + a*(-1 + 2*pow(mu,2)))))*excited_state + 
   ((2*sqrt(1 - pow(a,2))*gm*pow(mu,3)*sqrt(1 - pow(mu,2)) - 6*sqrt(1 - pow(a,2))*gm*pow(mu,5)*sqrt(1 - pow(mu,2)) + 
     4*sqrt(1 - pow(a,2))*gm*pow(mu,7)*sqrt(1 - pow(mu,2)) + 4*pow(a,6)*sqrt(1 - pow(a,2))*gm*mu*(1 - 2*pow(mu,2))*sqrt(1 - pow(mu,2)) + 
     2*pow(a,8)*sqrt(1 - pow(a,2))*gm*mu*sqrt(1 - pow(mu,2))*(-1 + 2*pow(mu,2)) + pow(a,9)*gm*(1 + 4*pow(mu,2) - 4*pow(mu,4)) - 
     4*pow(a,2)*sqrt(1 - pow(a,2))*gm*pow(mu,3)*sqrt(1 - pow(mu,2))*(1 - 3*pow(mu,2) + 2*pow(mu,4)) + pow(a,7)*gm*(-3 - 10*pow(mu,2) + 10*pow(mu,4)) + 
     3*a*gm*pow(mu,2)*(1 - 3*pow(mu,2) + 4*pow(mu,4) - 2*pow(mu,6)) + 
     2*pow(a,4)*sqrt(1 - pow(a,2))*gm*mu*sqrt(1 - pow(mu,2))*(-1 + 3*pow(mu,2) - 3*pow(mu,4) + 2*pow(mu,6)) + 
     pow(a,5)*gm*(3 + 11*pow(mu,2) - 15*pow(mu,4) + 8*pow(mu,6) - 4*pow(mu,8)) + 
     pow(a,3)*gm*(-1 - 8*pow(mu,2) + 18*pow(mu,4) - 20*pow(mu,6) + 10*pow(mu,8)) - 
     sqrt(pow(-1 + pow(a,2),4)*pow(gm,2)*(4*pow(a,9)*sqrt(1 - pow(a,2))*mu*(1 - 2*pow(mu,2))*sqrt(1 - pow(mu,2)) - 
         4*pow(mu,6)*pow(-1 + pow(mu,2),3) - 12*a*sqrt(1 - pow(a,2))*pow(mu,5)*pow(1 - pow(mu,2),2.5)*(-1 + 2*pow(mu,2)) + 
         pow(a,10)*(1 - 8*pow(mu,2) + 8*pow(mu,4)) + 9*pow(a,2)*pow(mu,4)*pow(-1 + pow(mu,2),2)*(1 - 8*pow(mu,2) + 8*pow(mu,4)) - 
         4*pow(a,7)*sqrt(1 - pow(a,2))*mu*sqrt(1 - pow(mu,2))*(1 + 6*pow(mu,2) - 24*pow(mu,4) + 16*pow(mu,6)) - 
         32*pow(a,5)*sqrt(1 - pow(a,2))*pow(mu,3)*sqrt(1 - pow(mu,2))*(-1 + pow(mu,2) + 6*pow(mu,4) - 10*pow(mu,6) + 4*pow(mu,8)) + 
         4*pow(a,3)*sqrt(1 - pow(a,2))*pow(mu,3)*sqrt(1 - pow(mu,2))*(-1 - 13*pow(mu,2) + 62*pow(mu,4) - 80*pow(mu,6) + 32*pow(mu,8)) + 
         2*pow(a,8)*(-1 + 2*pow(mu,2) + 30*pow(mu,4) - 64*pow(mu,6) + 32*pow(mu,8)) - 
         2*pow(a,4)*pow(mu,2)*(3 - 7*pow(mu,2) - 88*pow(mu,4) + 284*pow(mu,6) - 288*pow(mu,8) + 96*pow(mu,10)) + 
         pow(a,6)*(1 + 10*pow(mu,2) - 90*pow(mu,4) + 32*pow(mu,6) + 304*pow(mu,8) - 384*pow(mu,10) + 128*pow(mu,12)))))/
   (2.*pow(pow(a,2) - pow(a,4) - pow(mu,2) + pow(mu,4),2)))*ground_state;
}

Vector2cd Phi_phi (const VectorXcd &psi, double t) {
  double b = beta(t), gm = gamma_m(t), mu = real(psi(1));
  if (abs(mu) < .5)
    return sqrt(1.-mu*mu)*excited_state +
      ((2*b*(dt + dt*gm*(-1 + pow(mu,2))))/sqrt(1 - pow(mu,2)) + 
      mu*(1 + gm*(-2 + 3*dt + 2*pow(mu,2) - 2*dt*pow(mu,2)) + dt*pow(gm,2)*(-3 + 5*pow(mu,2) - 2*pow(mu,4))))*ground_state;
  return ((2*b*dt*(-1 + gm*(-1 + pow(mu,2))))/mu + sqrt(1 - pow(mu,2))*
    (1 + dt*pow(gm,2)*(-3 + 5*pow(mu,2) - 2*pow(mu,4)) + gm*(2 - 2*pow(mu,2) + dt*(-3 + 2*pow(mu,2)))))*excited_state +
    mu*ground_state;
}

Vector2cd Phi_phi_perp  (const VectorXcd &psi, double t) {
  double b = beta(t), gm = gamma_m(t), mu = real(psi(1));
  if (abs(mu) < .5)
    return (-((gm*pow(mu,2)*((2 + 3*dt*gm)*mu - (2 + 5*dt*gm)*pow(mu,3) + 2*dt*gm*pow(mu,5) + 2*b*dt*sqrt(1 - pow(mu,2))))/(-1 + pow(mu,2))))*excited_state;
  return mu*excited_state + 
    (b*dt*(2/mu - 2*gm*mu) + sqrt(1 - pow(mu,2))*(-1 + 2*gm*pow(mu,2) + dt*gm*(-3 + 2*pow(mu,2))*(-1 + gm*pow(mu,2))))*ground_state;
}

double observable (const MatrixXcd &rho) {return real((rho*sigma_z).trace());}

int main () {
  int N_ensemble = 1000, Ncopies = 1, dimH = 2, Ntraj = 10, Npoints = (int)((tmax-tmin)/dt);
  bool printTraj = true;

  Vector2cd initialState;
  initialState << 0.8, 0.7;

  // First thing: generating the deterministic trajectories
  int i = 1;
  phi[0] = ground_state;
  phi_perp[0] = excited_state;
  psi[0] = initialState.normalized();
  for (double t = tmin; t < tmax; t += dt) {
    MatrixXcd K = H(t) - .5*I*Gamma(t);

    phi[i] = phi[i-1] - I*dt*K*phi[i-1] - .5*dt*Phi_phi(phi[i-1],t);
    phi[i].normalize();
    if (real(phi[i](0)) < 0. || (real(phi[i](0)) < .01 && real(phi[i](1)) < 0.)) phi[i] *= -1.;

    phi_perp[i] = phi_perp[i-1] - I*dt*K*phi_perp[i-1] - .5*dt*Phi_phi_perp(phi[i-1],t);
    phi_perp[i].normalize();
    if (real(phi_perp[i](0)) < 0. || (real(phi_perp[i](0)) < .01 && real(phi_perp[i](1)) < 0.)) phi_perp[i] *= -1.;

    psi[i] = psi[i-1] - I*dt*K*psi[i-1] - .5*dt*Phi_psi(psi[i-1],t,real(phi[i](1)));
    psi[i].normalize();
    if (real(psi[i](0)) < 0. || (real(psi[i](0)) < .01 && real(psi[i](1)) < 0.)) psi[i] *= -1.;

    i++;
  }

  ofstream out;
  out.open("driven_gamma_minus_ens.txt");
  out << Npoints << endl;
  out << tmax << endl;
  out << dt << endl;

  // Now the jumps (plus exact solution)
  int N = 10000, N_psi = N, N_phi = 0, N_phi_perp = 0;
  int N_psi_old, N_phi_old, N_phi_perp_old;
  i = 0;

  Matrix2cd rho_ex = projector(psi[0]);

  for (double t = tmin; t < tmax; t += dt) {
    // Printing the observables
    Matrix2cd rho = ((double)N_psi/(double)N)*projector(psi[i]) + ((double)N_phi/(double)N)*projector(phi[i]) + ((double)N_phi_perp/(double)N)*projector(phi_perp[i]);
    out << observable(rho) << " " << observable(projector(psi[i])) << " " << observable(projector(phi[i])) << " " << observable(projector(phi_perp[i])) << " " << observable(rho_ex) << endl;

    rho_ex += (-I*comm(H(t),rho_ex) + J(rho_ex,t) - .5*anticomm(Gamma(t),rho_ex))*dt;

    N_psi_old = N_psi;
    N_phi_old = N_phi;
    N_phi_perp_old = N_phi_perp;

    // In psi
    Matrix2cd R = J(projector(psi[i]),t) + .5*(psi[i]*Phi_psi(psi[i],t,real(phi[i+1](1))).adjoint() + Phi_psi(psi[i],t,real(phi[i+1](1)))*psi[i].adjoint());
    ComplexEigenSolver<MatrixXcd> eigs;
    eigs.compute(R);
    VectorXcd eigval = eigs.eigenvalues();
    MatrixXcd eigvec = eigs.eigenvectors();
    double l_phi, l_perp;

    if ((eigval.col(0).normalized() - phi[i+1]).norm() < .1) {
      l_phi = real(eigval(0));
      l_perp = real(eigval(1));
    } else {
      l_perp = real(eigval(0));
      l_phi = real(eigval(1));
    }
    if (l_phi < 0. || l_perp < 0.)
      cout << "phi, t = " << t << ", mu = " << real(phi[i](1)) << "\teigs: " << l_phi << ", " << l_perp << endl;
    l_phi = norm(l_phi)*dt;
    l_perp= norm(l_perp)*dt;

    for (int j = 0; j < N_psi_old; ++j) {
      double z = (double)rand()/RAND_MAX;
      if (z <= l_phi) {
        N_psi--;
        N_phi++;
      } else if (z <= l_phi+l_perp) {
        N_psi--;
        N_phi_perp++;
      }
    }


    // In phi
    R = J(projector(phi[i]),t) + .5*(phi[i]*Phi_phi(phi[i],t).adjoint() + Phi_phi(phi[i],t)*phi[i].adjoint());
    eigs.compute(R);
    eigval = eigs.eigenvalues();
    eigvec = eigs.eigenvectors();

    if ((eigval.col(0).normalized() - phi[i+1]).norm() < .1) {
      l_phi = real(eigval(0));
      l_perp = real(eigval(1));
    } else {
      l_perp = real(eigval(0));
      l_phi = real(eigval(1));
    }
    if (l_phi < 0. || l_perp < 0.)
      cout << "phi, t = " << t << ", mu = " << real(phi[i](1)) << "\teigs: " << l_phi << ", " << l_perp << endl;
    l_phi = norm(l_phi)*dt;
    l_perp= norm(l_perp)*dt;

    for (int j = 0; j < N_phi_old; ++j) {
      double z = (double)rand()/RAND_MAX;
      if (z <= l_perp) {
        N_phi--;
        N_phi_perp++;
      }
    }

    // In phi_perp
    R = J(projector(phi_perp[i]),t) + .5*(phi_perp[i]*Phi_phi_perp(phi[i],t).adjoint() + Phi_phi_perp(phi[i],t)*phi_perp[i].adjoint());
    eigs.compute(R);
    eigval = eigs.eigenvalues();
    eigvec = eigs.eigenvectors();

    if ((eigval.col(0).normalized() - phi[i+1]).norm() < .1) {
      l_phi = real(eigval(0));
      l_perp = real(eigval(1));
    } else {
      l_perp = real(eigval(0));
      l_phi = real(eigval(1));
    }
    if (l_phi < 0. || l_perp < 0.)
      cout << "phi_perp, t = " << t << ", mu = " << real(phi[i](1)) << "\teigs: " << l_phi << ", " << l_perp << endl;
    l_phi = norm(l_phi)*dt;
    l_perp= norm(l_perp)*dt;

    for (int j = 0; j < N_phi_perp_old; ++j) {
      double z = (double)rand()/RAND_MAX;
      if (z <= l_phi) {
        N_phi_perp--;
        N_phi++;
      }
    }

    ++i;
  }

  qubit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj, true, threshold);
  jump.set_initial_state(initialState);

  return 0;
}