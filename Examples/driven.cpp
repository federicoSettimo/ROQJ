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

double gamma_m (double t) {return 1.;}
double gamma_p (double t) {return .9;}
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

Matrix2cd comm (const Matrix2cd &A, const Matrix2cd &B) {return A*B - B*A;}

Matrix2cd anticomm (const Matrix2cd &A, const Matrix2cd &B) {return A*B + B*A;}

int main () {
  int Nensemble = 10000, Npsi = Nensemble, Ng = 0, Ne = 0, Np = 0, Nm = 0;
  double tmin = 0., tmax = 10., dt = 0.01, threshold = 1e-3;
  double pe, pg, pp, pm;

  ofstream out;
  out.open("driven.txt");
  out << (int)(tmax/dt) << endl << tmax << endl << dt << endl;

  Vector2cd psi;
  psi << .7, .8;
  psi.normalize();

  Matrix2cd exact = projector(psi), rho;

  for (double t = 0.; t < tmax; t += dt) {
    // Old populations
    int Npsi_old = Npsi, Ng_old = Ng, Ne_old = Ne, Np_old = Np, Nm_old = Nm;
    // Rates and parameters
    double a = real(psi(1)), gp = gamma_p(t), gm = gamma_m(t), a2 = a*a, b = beta(t);
    // Fractions of states
    double fpsi = (double)Npsi/(double)Nensemble, fp = (double)Np/(double)Nensemble, fm = (double)Nm/(double)Nensemble, fg = (double)Ng/(double)Nensemble, fe = (double)Ne/(double)Nensemble;
    rho = fpsi*projector(psi) + fg*projector(ground_state) + fe*projector(excited_state) + fp*projector(plus_state) + fm*projector(minus_state);
    out << observable(rho) << " " << observable(exact) << " " << observable(projector(psi)) << endl;
    out << fpsi << " " << fg << " " << fe << " " << fp << " " << fm << endl;
    out << gp << " " << gm << " " << b << endl;

    // Updating exact
    exact += dt*(-I*comm(H(t),exact) + J(exact,t) - .5*anticomm(Gamma(t),exact));

    // Direct jumps from psi
    pe = a2*gp*dt;
    pg = (1.-a2)*gm*dt;
    for (int i = 0; i < Npsi_old; ++i) {
      double z = ((double)rand()/(double)RAND_MAX);
      if (z < pe) {Npsi--; Ne++;}
      else if (z < pe+pg) {Npsi--; Ng++;}
    }
    // Update psi
    Matrix2cd K = H(t) - .5*I*Gamma(t);
    psi -= I*dt*K*psi;
    psi.normalize();
    if (pe < 0. || pg < 0.) cerr << "Error in psi, t = " << t << "\t" << pe << ", " << pg << endl;

    // From g (positive rates for b < gp/2)
    pp = -2.*b + gp;
    pm = 2.*b + gp;
    for (int i = 0; i < Ng_old; ++i) {
      double z = ((double)rand()/(double)RAND_MAX);
      if (z < pp) {Ng--; Np++;}
      else if (z < pp+pm) {Ng--; Np++;}
    }
    if (pp < 0. || pm < 0.) cerr << "Error in g, t = " << t << "\t" << pp << ", " << pm << endl;

    // From e (positive rates for b < gm/2)
    pp = 2.*b + gm;
    pm = -2.*b + gm;
    for (int i = 0; i < Ne_old; ++i) {
      double z = ((double)rand()/(double)RAND_MAX);
      if (z < pp) {Ne--; Np++;}
      else if (z < pp+pm) {Ne--; Np++;}
    }
    if (pp < 0. || pm < 0.) cerr << "Error in e, t = " << t << "\t" << pp << ", " << pm << endl;

    // From + (positive rates for 2b < gp - gm/2)
    pg = 2.*b + gm - .5*gp;
    pe = -2.*b + gp - .5*gm;
    for (int i = 0; i < Np_old; ++i) {
      double z = ((double)rand()/(double)RAND_MAX);
      if (z < pg) {Np--; Ng++;}
      else if (z < pg+pe) {Np--; Ne++;}
    }
    if (pg < 0. || pe < 0.) cerr << "Error in +, t = " << t << "\t" << pg << ", " << pe << endl;

    // From - (positive rates for 2b < gp - gm/2)
    pg = -2.*b + gm - .5*gp;
    pe = 2.*b + gp - .5*gm;
    for (int i = 0; i < Nm_old; ++i) {
      double z = ((double)rand()/(double)RAND_MAX);
      if (z < pg) {Nm--; Ng++;}
      else if (z < pg+pe) {Nm--; Ne++;}
    }
    if (pg < 0. || pe < 0.) cerr << "Error in -, t = " << t << "\t" << pg << ", " << pe << endl;
  }

  return 0;
}