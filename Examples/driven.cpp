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
/*double gamma_m (double t) {return 1.;}
double gamma_p (double t) {return .9;}
double beta (double t) {
  double gp = gamma_p(t), gm = gamma_m(t);
  return .5*min(gp,gm);
}*/
double gamma_m (double t) {return 1.;}
//double gamma_p (double t) {return .5*(tanh(.5*t)+1);}
double gamma_p (double t) {return 0.;}
double beta (double t) {return .2;}

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

int Nensemble = 10000;
int Npsi = Nensemble, Ng = 0, Ne = 0, Np = 0, Nm = 0;
int Npsi_old, Ng_old, Ne_old, Np_old, Nm_old;
double tmin = 0., tmax = 10., dt = 0.001, threshold = 1e-3;

/*void jump (double p1, double p2, int &N, int &N1, int &N2, int Nold, int N1old, int N2old) {
  double z = ((double)rand()/(double)RAND_MAX);
  if (p1 >= 0. && p2 >= 0.) { // Direct x -> 1,2
    if (z < p1) {N--; N1++;}
    else if (z < p1+p2) {N--; N2++;}
  }
  else if (p1 >= 0. && p2 < 0.) { // Direct x -> 1, reverse 2 -> x
    if (z < p1) {N--; N1++;}
    else if (N2 != 0) {
      double p_rev_2 = abs(p2)/(double)N2old*(double)Nold;
      if (z < p1 + p_rev_2) {N++; N2--;}
    }
  }
  else if (p2 >= 0. && p1 < 0.) { // Direct x -> 2, reverse 1 -> x
    if (z < p2) {N--; N2++;}
    else if (N1 != 0) {
      double p_rev_1 = abs(p1)/(double)N1old*(double)Nold;
      if (z < p2 + p_rev_1) {N++; N1--;}
    }
  }
  else { // Reverse 1,2 -> x
    if (N1 != 0) {
      double p_rev_1 = -p1/(double)N1old*(double)Nold;
      if (z < p_rev_1) {N++; N1--;}
    }
    if (N2 != 0) {
      double p_rev_2 = -p2/(double)N2old*(double)Nold;
      if (z > 1. - p_rev_2) {N++, N1--;}
    }
  }
}*/

void jump (double p1, double p2, int &N, int &N1, int &N2, int Nold, int N1old, int N2old) {
  double z = ((double)rand()/(double)RAND_MAX);
  if (p1 >= 0. && p2 >= 0.) { // Direct x -> 1,2
    if (z < p1) {N--; N1++;}
    else if (z < p1+p2) {N--; N2++;}
  }
  else if (p1 >= 0. && p2 < 0.) { // Direct x -> 1, reverse 2 -> x
    double p_rev_2 = -p2*(double)N2old/(double)Nold;
    if (z < p1) {N--; N1++;}
    else if (z < p1 + p_rev_2) {N++; N2--;}
  }
  else if (p2 >= 0. && p1 < 0.) { // Direct x -> 2, reverse 1 -> x
    double p_rev_1 = -p1*(double)N1old/(double)Nold;
    if (z < p2) {N--; N2++;}
    else if (z < p2 + p_rev_1) {N++; N1--;}
  }
  else { // Reverse 1,2 -> x
    double p_rev_1 = -p1*(double)N1old/(double)Nold, p_rev_2 = -p2*(double)N2old/(double)Nold;
    if (z < p_rev_1) {N++; N1--;}
    else if (z < p_rev_1 + p_rev_2) {N++, N1--;}
  }
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

  Vector2cd psi;
  psi << .7, .8;
  psi.normalize();

  Matrix2cd exact = projector(psi), rho;

  for (double t = 0.; t < tmax; t += dt) {
    // Old populations
    Npsi_old = Npsi; Ng_old = Ng; Ne_old = Ne; Np_old = Np; Nm_old = Nm;
    // Rates and parameters
    double a = real(psi(1)), gp = gamma_p(t), gm = gamma_m(t), a2 = a*a, b = beta(t);
    // Fractions of states
    double fpsi = (double)Npsi/(double)Nensemble, fp = (double)Np/(double)Nensemble, fm = (double)Nm/(double)Nensemble, fg = (double)Ng/(double)Nensemble, fe = (double)Ne/(double)Nensemble;

    rho = fpsi*projector(psi) + fg*projector(ground_state) + fe*projector(excited_state) + fp*projector(plus_state) + fm*projector(minus_state);
    out << observable(rho) << " " << observable(exact) << " " << observable(projector(psi)) << endl;
    out << coherence(rho) << " " << coherence(exact) << " " << coherence(projector(psi)) << endl;
    out << fpsi << " " << fg << " " << fe << " " << fp << " " << fm << endl;
    out << gp << " " << gm << " " << b << endl;

    // Updating exact
    //exact += dt*(-I*comm(H(t),exact) + J(exact,t) - .5*anticomm(Gamma(t),exact));
    // Updating the exact (with Runge-Kutta 4)
    exact += RK4(exact,t);

    // Updating psi
    Matrix2cd K = H(t) - .5*I*Gamma(t);
    psi -= I*dt*K*psi;
    psi.normalize();

    double psi_e = a2*gp*dt, psi_g = (1.-a2)*gm*dt;
    double g_p = (-2.*b + gp)*dt, g_m = (2.*b + gp)*dt; 
    double e_p = (2.*b + gm)*dt, e_m = (-2.*b + gm)*dt;
    double p_g = (2.*b + gm - .5*gp)*dt, p_e = (-2.*b + gp - .5*gm)*dt;
    double m_g = (-2.*b + gm - .5*gp)*dt, m_e = (2.*b + gp - .5*gm)*dt;
    out << psi_e << " " << psi_g << " " << g_p << " " << g_m << " " << e_p << " " << e_m << " " << p_g << " " << p_e << " " << m_g << " " << m_e << endl;

    for (int i = 0; i < Npsi_old; ++i) // psi
      jump(psi_e, psi_g, Npsi, Ne, Ng, Npsi_old, Ne_old, Ng_old);
    for (int i = 0; i < Ng_old; ++i) // g
      jump(g_p, g_m, Ng, Np, Nm, Ng_old, Np_old, Nm_old);
    for (int i = 0; i < Ne_old; ++i) // e
      jump(e_p, e_m, Ne, Np, Nm, Ne_old, Np_old, Nm_old);
    for (int i = 0; i < Np_old; ++i) // +
      jump(p_e, p_g, Np, Ne, Ng, Np_old, Ne_old, Ng_old);
    for (int i = 0; i < Nm_old; ++i) // -
      jump(m_e, m_g, Nm, Ne, Ng, Nm_old, Ne_old, Ng_old);
  }

  return 0;
}