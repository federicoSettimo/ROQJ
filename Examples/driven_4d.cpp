// Driven with 4 ensemble: psi, psi_perp, (psi + psi_perp)/sqrt(2), (psi - psi_perp)/sqrt(2)
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
double gamma_p (double t) {return .5*(sin(t)+1.);}
double beta (double t) {
  double gp = gamma_p(t), gm = gamma_m(t);
  return t > 1. ? .5*min(gp,gm) : 0.;
  return .5*min(gp,gm);
}

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
int Ncopies = 10, Nensemble = 1000, Npsi = Nensemble, Nperp = 0, Nplus = 0, Nmin = 0;
int Npsi_old, Nperp_old, Nplus_old, Nmin_old;

void jump (double p1, double p2, int &N, int &N1, int &N2) {
  double z = (double)rand()/(double)RAND_MAX;
  if (z < p1) {N--; N1++;}
  else if (z < p1+p2) {N--; N2++;}
  return;
}

double rev_jump (double l, int N, int N_target) {
  if (l > 0. || N == 0.) return 0.;
  return l*(double)N_target/(double)N;
}

Matrix2cd L (const Matrix2cd &rho, double t) {return -I*comm(H(t),rho) + J(rho,t) - .5*anticomm(Gamma(t),rho);}

// Runge-Kutta 4
Matrix2cd RK4 (const Matrix2cd &rho, double t) {
  Matrix2cd k1 = dt * L(rho, t), k2 = dt * L(rho + .5*k1, t + .5*dt), k3 = dt * L(rho + .5*k2, t + .5*dt), k4 = dt * L(rho + k3, t + dt);
  return (1./6.) * (k1 + 2.*k2 + 2.*k3 + k4);
}

Vector2cd Phi (const Vector2cd &psi, double t) {
  double mu = real(psi(1)), gp = gamma_p(t), gm = gamma_m(t), b = beta(t);
  return (-(((gp*pow(mu,2) + gm*(-1 + pow(mu,2)))*(gp*pow(mu,3)*pow(1 - 2*pow(mu,2),2)*
          (1 + (1 - 8*b*dt + 5*dt*gp)*pow(mu,2) + 2*(-7 + 4*b*dt - 4*dt*gp)*pow(mu,4) + 2*(10 + dt*gp)*pow(mu,6) - 8*(1 + dt*gp)*pow(mu,8) + 
            8*dt*gp*pow(mu,10) + (4 - 2*b*dt + dt*gp)*mu*sqrt(1 - pow(mu,2)) - 2*(6 + 4*b*dt - 3*dt*gp)*pow(mu,3)*sqrt(1 - pow(mu,2)) + 
            4*(2 + 2*b*dt - 3*dt*gp)*pow(mu,5)*sqrt(1 - pow(mu,2))) - 
         dt*pow(gm,2)*mu*(-1 + pow(mu,2))*(2 + 48*pow(mu,2) - 164*pow(mu,4) + 106*pow(mu,6) + 168*pow(mu,8) - 288*pow(mu,10) + 160*pow(mu,12) - 
            32*pow(mu,14) + 17*mu*sqrt(1 - pow(mu,2)) + 26*pow(mu,3)*sqrt(1 - pow(mu,2)) - 244*pow(mu,5)*sqrt(1 - pow(mu,2)) + 
            392*pow(mu,7)*sqrt(1 - pow(mu,2)) - 256*pow(mu,9)*sqrt(1 - pow(mu,2)) + 64*pow(mu,11)*sqrt(1 - pow(mu,2))) - 
         gm*((13 - 4*b*dt + 2*dt*gp)*mu + (-48 - 48*b*dt + 33*dt*gp)*pow(mu,3) + (1 + 188*b*dt - 145*dt*gp)*pow(mu,5) + (234 - 248*b*dt + 202*dt*gp)*pow(mu,7) + 
            4*(-114 + 36*b*dt - 41*dt*gp)*pow(mu,9) - 8*(-50 + 4*b*dt - 31*dt*gp)*pow(mu,11) - 16*(11 + 23*dt*gp)*pow(mu,13) + 32*(1 + 8*dt*gp)*pow(mu,15) - 
            64*dt*gp*pow(mu,17) + 2*sqrt(1 - pow(mu,2)) + (20 - 26*b*dt + 15*dt*gp)*pow(mu,2)*sqrt(1 - pow(mu,2)) + 
            2*(-71 + 9*b*dt - 4*dt*gp)*pow(mu,4)*sqrt(1 - pow(mu,2)) + 4*(78 + 26*b*dt - 25*dt*gp)*pow(mu,6)*sqrt(1 - pow(mu,2)) - 
            4*(80 + 48*b*dt - 35*dt*gp)*pow(mu,8)*sqrt(1 - pow(mu,2)) + 32*(5 + 4*b*dt + dt*gp)*pow(mu,10)*sqrt(1 - pow(mu,2)) - 
            16*(2 + 2*b*dt + 9*dt*gp)*pow(mu,12)*sqrt(1 - pow(mu,2)) + 64*dt*gp*pow(mu,14)*sqrt(1 - pow(mu,2)))))/
     (pow(1 - pow(mu,2),1.5)*pow(mu + sqrt(1 - pow(mu,2)),4)*
       (gp*pow(mu,3)*pow(1 - 2*pow(mu,2),2) + gm*(-1 + pow(mu,2))*(5*mu - 8*pow(mu,3) + 4*pow(mu,5) + 2*sqrt(1 - pow(mu,2)))))))*excited_state;
}

Vector2cd PhiPerp (const Vector2cd &psi, double t) {
  double mu = real(psi(1)), gp = gamma_p(t), gm = gamma_m(t), b = beta(t);
  return (-((2*b*dt*(-1 + pow(mu,2))*(gp*(4 - 5*pow(mu,2)) + 3*gm*(-1 + pow(mu,2)))*
        (gm*mu*(3 + 33*pow(mu,2) - 60*pow(mu,4) + 40*pow(mu,6) - 32*pow(mu,8) + 16*pow(mu,10) + 18*mu*sqrt(1 - pow(mu,2)) + 
             8*pow(mu,3)*sqrt(1 - pow(mu,2)) - 24*pow(mu,5)*sqrt(1 - pow(mu,2))) + 
          gp*(8*mu - 63*pow(mu,3) + 80*pow(mu,5) + 24*pow(mu,7) - 64*pow(mu,9) + 16*pow(mu,11) + 2*sqrt(1 - pow(mu,2)) - 
             4*pow(mu,2)*sqrt(1 - pow(mu,2)) - 72*pow(mu,4)*sqrt(1 - pow(mu,2)) + 144*pow(mu,6)*sqrt(1 - pow(mu,2)) - 64*pow(mu,8)*sqrt(1 - pow(mu,2))))
         - (gm - gp)*mu*(3 - 7*pow(mu,2) + 4*pow(mu,4))*(gm*mu*(-1 + pow(mu,2))*
           (-18*mu - 8*pow(mu,3) + 24*pow(mu,5) - 3*sqrt(1 - pow(mu,2)) - 36*pow(mu,2)*sqrt(1 - pow(mu,2)) + 24*pow(mu,4)*sqrt(1 - pow(mu,2)) - 
             16*pow(mu,6)*sqrt(1 - pow(mu,2)) + 16*pow(mu,8)*sqrt(1 - pow(mu,2))) + 
          gp*(2 - 6*pow(mu,2) - 68*pow(mu,4) + 216*pow(mu,6) - 208*pow(mu,8) + 64*pow(mu,10) + 8*mu*sqrt(1 - pow(mu,2)) - 
             63*pow(mu,3)*sqrt(1 - pow(mu,2)) + 80*pow(mu,5)*sqrt(1 - pow(mu,2)) + 24*pow(mu,7)*sqrt(1 - pow(mu,2)) - 
             64*pow(mu,9)*sqrt(1 - pow(mu,2)) + 16*pow(mu,11)*sqrt(1 - pow(mu,2)))) - 
       dt*(gp*(4 - 5*pow(mu,2)) + 3*gm*(-1 + pow(mu,2)))*(pow(gm,2)*mu*(-1 + pow(mu,2))*
           (-3 - 63*pow(mu,2) + 182*pow(mu,4) - 116*pow(mu,6) - 32*pow(mu,10) + 32*pow(mu,12) - 24*mu*sqrt(1 - pow(mu,2)) - 
             32*pow(mu,3)*sqrt(1 - pow(mu,2)) + 226*pow(mu,5)*sqrt(1 - pow(mu,2)) - 248*pow(mu,7)*sqrt(1 - pow(mu,2)) + 
             144*pow(mu,9)*sqrt(1 - pow(mu,2)) - 96*pow(mu,11)*sqrt(1 - pow(mu,2)) + 32*pow(mu,13)*sqrt(1 - pow(mu,2))) + 
          gm*gp*(-1 + pow(mu,2))*(-9*mu + 140*pow(mu,3) - 200*pow(mu,5) - 332*pow(mu,7) + 832*pow(mu,9) - 560*pow(mu,11) + 128*pow(mu,13) - 
             2*sqrt(1 - pow(mu,2)) + 13*pow(mu,2)*sqrt(1 - pow(mu,2)) + 214*pow(mu,4)*sqrt(1 - pow(mu,2)) - 644*pow(mu,6)*sqrt(1 - pow(mu,2)) + 
             456*pow(mu,8)*sqrt(1 - pow(mu,2)) + 80*pow(mu,10)*sqrt(1 - pow(mu,2)) - 192*pow(mu,12)*sqrt(1 - pow(mu,2)) + 
             64*pow(mu,14)*sqrt(1 - pow(mu,2))) + pow(gp,2)*(-10*mu + 93*pow(mu,3) - 213*pow(mu,5) + 114*pow(mu,7) + 48*pow(mu,9) + 112*pow(mu,11) - 
             240*pow(mu,13) + 96*pow(mu,15) - 2*sqrt(1 - pow(mu,2)) + 2*pow(mu,2)*sqrt(1 - pow(mu,2)) + 119*pow(mu,4)*sqrt(1 - pow(mu,2)) - 
             416*pow(mu,6)*sqrt(1 - pow(mu,2)) + 490*pow(mu,8)*sqrt(1 - pow(mu,2)) - 256*pow(mu,10)*sqrt(1 - pow(mu,2)) + 
             160*pow(mu,12)*sqrt(1 - pow(mu,2)) - 128*pow(mu,14)*sqrt(1 - pow(mu,2)) + 32*pow(mu,16)*sqrt(1 - pow(mu,2)))))/
     (pow(1 - pow(mu,2),1.5)*pow(mu + sqrt(1 - pow(mu,2)),4)*(1 + 2*mu*sqrt(1 - pow(mu,2)))*
       (gm*mu*(3 - 3*pow(mu,2) + 2*pow(mu,3)*sqrt(1 - pow(mu,2))) + 
         gp*(-4*mu + 5*pow(mu,3) + 2*sqrt(1 - pow(mu,2)) - 4*pow(mu,2)*sqrt(1 - pow(mu,2)) + 2*pow(mu,4)*sqrt(1 - pow(mu,2)))))))*excited_state +
    ((gp + 2*gm*pow(mu,2) - 6*gp*pow(mu,2) - 2*gm*pow(mu,4) + 6*gp*pow(mu,4))/sqrt(1 - pow(mu,2)))*ground_state;
}

Vector2cd PhiPlus (const Vector2cd &psi, double t) {
  double mu = real(psi(1)), gp = gamma_p(t), gm = gamma_m(t), b = beta(t);
  return ((gm*(8*mu + 44*pow(mu,3) - 196*pow(mu,5) + 224*pow(mu,7) - 80*pow(mu,9) + sqrt(1 - pow(mu,2)) + 30*pow(mu,2)*sqrt(1 - pow(mu,2)) - 
        24*pow(mu,4)*sqrt(1 - pow(mu,2)) - 120*pow(mu,6)*sqrt(1 - pow(mu,2)) + 176*pow(mu,8)*sqrt(1 - pow(mu,2)) - 64*pow(mu,10)*sqrt(1 - pow(mu,2))) - 
     gp*(6*mu - 4*pow(mu,3) + 44*pow(mu,5) - 128*pow(mu,7) + 80*pow(mu,9) + sqrt(1 - pow(mu,2)) + 10*pow(mu,2)*sqrt(1 - pow(mu,2)) + 
        24*pow(mu,6)*sqrt(1 - pow(mu,2)) - 112*pow(mu,8)*sqrt(1 - pow(mu,2)) + 64*pow(mu,10)*sqrt(1 - pow(mu,2))))/
   (2.*sqrt(2 - 2*pow(mu,2))*pow(mu + sqrt(1 - pow(mu,2)),3)*(1 + 2*mu*sqrt(1 - pow(mu,2)))))*excited_state + 
    ((8*dt*pow(gm,3)*pow(mu,2)*pow(-1 + pow(mu,2),4)*(-3 + 4*pow(mu,2))*
      (13*mu + 244*pow(mu,3) + 56*pow(mu,5) - 2184*pow(mu,7) + 2896*pow(mu,9) - 384*pow(mu,11) - 1024*pow(mu,13) + 384*pow(mu,15) + 
        sqrt(1 - pow(mu,2)) + 76*pow(mu,2)*sqrt(1 - pow(mu,2)) + 448*pow(mu,4)*sqrt(1 - pow(mu,2)) - 984*pow(mu,6)*sqrt(1 - pow(mu,2)) - 
        944*pow(mu,8)*sqrt(1 - pow(mu,2)) + 2688*pow(mu,10)*sqrt(1 - pow(mu,2)) - 1408*pow(mu,12)*sqrt(1 - pow(mu,2)) + 
        128*pow(mu,14)*sqrt(1 - pow(mu,2))) + pow(gm,2)*pow(-1 + pow(mu,2),2)*
      ((19 - 24*b*dt)*mu - 2*(-311 + 884*b*dt + 102*dt*gp)*pow(mu,3) - 8*(15 + 816*b*dt + 503*dt*gp)*pow(mu,5) + 8*(-2033 + 5784*b*dt + 2585*dt*gp)*pow(mu,7) - 
        32*(-1459 + 1462*b*dt - 370*dt*gp)*pow(mu,9) - 32*(707 + 2684*b*dt + 8278*dt*gp)*pow(mu,11) + 64*(-1320 + 3352*b*dt + 10579*dt*gp)*pow(mu,13) - 
        128*(-1189 + 1312*b*dt + 5693*dt*gp)*pow(mu,15) + 256*(-401 + 204*b*dt + 950*dt*gp)*pow(mu,17) - 1024*(-28 + 4*b*dt - 167*dt*gp)*pow(mu,19) - 
        1024*(2 + 159*dt*gp)*pow(mu,21) + 36864*dt*gp*pow(mu,23) + sqrt(1 - pow(mu,2)) - 6*(-25 + 52*b*dt + 2*dt*gp)*pow(mu,2)*sqrt(1 - pow(mu,2)) - 
        8*(680*b*dt + 169*(-1 + dt*gp))*pow(mu,4)*sqrt(1 - pow(mu,2)) + 8*(-851 + 808*b*dt - 363*dt*gp)*pow(mu,6)*sqrt(1 - pow(mu,2)) + 
        32*(-173 + 1694*b*dt + 1700*dt*gp)*pow(mu,8)*sqrt(1 - pow(mu,2)) - 32*(-2015 + 4356*b*dt + 4740*dt*gp)*pow(mu,10)*sqrt(1 - pow(mu,2)) + 
        64*(-1736 + 1592*b*dt + 1101*dt*gp)*pow(mu,12)*sqrt(1 - pow(mu,2)) + 128*(503 + 96*b*dt + 2735*dt*gp)*pow(mu,14)*sqrt(1 - pow(mu,2)) - 
        256*(-41 + 164*b*dt + 2666*dt*gp)*pow(mu,16)*sqrt(1 - pow(mu,2)) + 1024*(-23 + 12*b*dt + 496*dt*gp)*pow(mu,18)*sqrt(1 - pow(mu,2)) - 
        3072*(-2 + 51*dt*gp)*pow(mu,20)*sqrt(1 - pow(mu,2)) + 12288*dt*gp*pow(mu,22)*sqrt(1 - pow(mu,2))) + 
     pow(gp,2)*(-19*mu + 8*(-95 + 3*b*dt)*pow(mu,3) + (-2971 + 2272*b*dt - 180*dt*gp)*pow(mu,5) + 2*(7663 + 4404*b*dt - 1854*dt*gp)*pow(mu,7) + 
        (9040 - 68640*b*dt + 6200*dt*gp)*pow(mu,9) + 8*(-13501 + 12688*b*dt + 5175*dt*gp)*pow(mu,11) + 64*(2365 + 963*b*dt - 1943*dt*gp)*pow(mu,13) - 
        32*(-353 + 10404*b*dt - 1858*dt*gp)*pow(mu,15) + 64*(-3638 + 6208*b*dt + 2841*dt*gp)*pow(mu,17) - 128*(-2049 + 1736*b*dt + 2211*dt*gp)*pow(mu,19) + 
        256*(-521 + 220*b*dt + 446*dt*gp)*pow(mu,21) - 1024*(-30 + 4*b*dt - 49*dt*gp)*pow(mu,23) - 1024*(2 + 53*dt*gp)*pow(mu,25) + 12288*dt*gp*pow(mu,27) - 
        sqrt(1 - pow(mu,2)) - 160*pow(mu,2)*sqrt(1 - pow(mu,2)) + (-2221 + 360*b*dt - 12*dt*gp)*pow(mu,4)*sqrt(1 - pow(mu,2)) + 
        2*(219 + 3708*b*dt - 574*dt*gp)*pow(mu,6)*sqrt(1 - pow(mu,2)) - 8*(-3438 + 1460*b*dt + 691*dt*gp)*pow(mu,8)*sqrt(1 - pow(mu,2)) - 
        8*(5103 + 8496*b*dt - 3883*dt*gp)*pow(mu,10)*sqrt(1 - pow(mu,2)) + 64*(-923 + 3521*b*dt - 170*dt*gp)*pow(mu,12)*sqrt(1 - pow(mu,2)) - 
        32*(-6571 + 7964*b*dt + 3448*dt*gp)*pow(mu,14)*sqrt(1 - pow(mu,2)) + 64*(-3302 + 1360*b*dt + 2463*dt*gp)*pow(mu,16)*sqrt(1 - pow(mu,2)) + 
        128*(515 + 440*b*dt + 161*dt*gp)*pow(mu,18)*sqrt(1 - pow(mu,2)) - 256*(-129 + 212*b*dt + 786*dt*gp)*pow(mu,20)*sqrt(1 - pow(mu,2)) + 
        1024*(-29 + 12*b*dt + 164*dt*gp)*pow(mu,22)*sqrt(1 - pow(mu,2)) - 3072*(-2 + 17*dt*gp)*pow(mu,24)*sqrt(1 - pow(mu,2)) + 
        4096*dt*gp*pow(mu,26)*sqrt(1 - pow(mu,2))) + 4*gm*gp*mu*(-1 + pow(mu,2))*
      (-1 - 15*(4 + 3*dt*gp)*pow(mu,2) + (635 - 978*dt*gp)*pow(mu,4) + (3142 + 466*dt*gp)*pow(mu,6) + 2*(-11086 + 7119*dt*gp)*pow(mu,8) - 
        8*(-4609 + 3144*dt*gp)*pow(mu,10) + (3184 - 39424*dt*gp)*pow(mu,12) + 144*(-554 + 1155*dt*gp)*pow(mu,14) - 32*(-3230 + 6019*dt*gp)*pow(mu,16) + 
        128*(-461 + 523*dt*gp)*pow(mu,18) + 256*(58 + 163*dt*gp)*pow(mu,20) - 256*(4 + 159*dt*gp)*pow(mu,22) + 9216*dt*gp*pow(mu,24) - 
        (13 + 3*dt*gp)*mu*sqrt(1 - pow(mu,2)) - (53 + 290*dt*gp)*pow(mu,3)*sqrt(1 - pow(mu,2)) + 2*(1303 - 863*dt*gp)*pow(mu,5)*sqrt(1 - pow(mu,2)) + 
        2*(-2226 + 3299*dt*gp)*pow(mu,7)*sqrt(1 - pow(mu,2)) + 8*(-2175 + 1156*dt*gp)*pow(mu,9)*sqrt(1 - pow(mu,2)) + 
        16*(4013 - 3371*dt*gp)*pow(mu,11)*sqrt(1 - pow(mu,2)) + 16*(-4850 + 2833*dt*gp)*pow(mu,13)*sqrt(1 - pow(mu,2)) + 
        32*(1002 + 2225*dt*gp)*pow(mu,15)*sqrt(1 - pow(mu,2)) - 128*(-85 + 1307*dt*gp)*pow(mu,17)*sqrt(1 - pow(mu,2)) + 
        1024*(-13 + 124*dt*gp)*pow(mu,19)*sqrt(1 - pow(mu,2)) - 768*(-4 + 51*dt*gp)*pow(mu,21)*sqrt(1 - pow(mu,2)) + 
        3072*dt*gp*pow(mu,23)*sqrt(1 - pow(mu,2)) - 2*b*dt*(-3 + 4*pow(mu,2))*
         (1 + 95*pow(mu,2) + 420*pow(mu,4) - 2572*pow(mu,6) + 2728*pow(mu,8) + 4256*pow(mu,10) - 11776*pow(mu,12) + 9792*pow(mu,14) - 3200*pow(mu,16) + 
           256*pow(mu,18) + 15*mu*sqrt(1 - pow(mu,2)) + 316*pow(mu,3)*sqrt(1 - pow(mu,2)) - 292*pow(mu,5)*sqrt(1 - pow(mu,2)) - 
           2952*pow(mu,7)*sqrt(1 - pow(mu,2)) + 7712*pow(mu,9)*sqrt(1 - pow(mu,2)) - 6144*pow(mu,11)*sqrt(1 - pow(mu,2)) - 
           320*pow(mu,13)*sqrt(1 - pow(mu,2)) + 2432*pow(mu,15)*sqrt(1 - pow(mu,2)) - 768*pow(mu,17)*sqrt(1 - pow(mu,2)))))/
   (2.*sqrt(2 - 2*pow(mu,2))*pow(mu + sqrt(1 - pow(mu,2)),2)*(2*mu - 2*pow(mu,3) + sqrt(1 - pow(mu,2)))*(1 + 2*mu*sqrt(1 - pow(mu,2)))*
     (gm*(-1 + pow(mu,2))*(1 + 28*pow(mu,2) - 12*pow(mu,4) - 144*pow(mu,6) + 192*pow(mu,8) - 64*pow(mu,10) + 8*mu*sqrt(1 - pow(mu,2)) + 
          48*pow(mu,3)*sqrt(1 - pow(mu,2)) - 120*pow(mu,5)*sqrt(1 - pow(mu,2)) + 32*pow(mu,7)*sqrt(1 - pow(mu,2)) + 32*pow(mu,9)*sqrt(1 - pow(mu,2))) + 
       gp*(-1 - 39*pow(mu,2) + 36*pow(mu,4) + 164*pow(mu,6) - 352*pow(mu,8) + 256*pow(mu,10) - 64*pow(mu,12) - 10*mu*sqrt(1 - pow(mu,2)) - 
          64*pow(mu,3)*sqrt(1 - pow(mu,2)) + 192*pow(mu,5)*sqrt(1 - pow(mu,2)) - 152*pow(mu,7)*sqrt(1 - pow(mu,2)) + 32*pow(mu,11)*sqrt(1 - pow(mu,2))))
     ))*ground_state;
}

Vector2cd PhiMin (const Vector2cd &psi, double t) {
  double mu = real(psi(1)), gp = gamma_p(t), gm = gamma_m(t), b = beta(t);
  return ((gm*(-4*mu + 8*pow(mu,3) - 4*pow(mu,5) + sqrt(1 - pow(mu,2)) + 10*pow(mu,2)*sqrt(1 - pow(mu,2)) - 28*pow(mu,4)*sqrt(1 - pow(mu,2)) + 
        16*pow(mu,6)*sqrt(1 - pow(mu,2))) + gp*(2*mu - 4*pow(mu,5) - sqrt(1 - pow(mu,2)) + 2*pow(mu,2)*sqrt(1 - pow(mu,2)) - 
        12*pow(mu,4)*sqrt(1 - pow(mu,2)) + 16*pow(mu,6)*sqrt(1 - pow(mu,2))))/(2.*sqrt(2 - 2*pow(mu,2))*(-mu + sqrt(1 - pow(mu,2)))))*excited_state +
    ((-8*dt*pow(gm,3)*pow(mu,2)*pow(-1 + pow(mu,2),4)*(-3 + 4*pow(mu,2))*
      (3*mu - 6*pow(mu,3) + 4*pow(mu,5) - sqrt(1 - pow(mu,2)) - 6*pow(mu,2)*sqrt(1 - pow(mu,2)) + 12*pow(mu,4)*sqrt(1 - pow(mu,2))) - 
     pow(gm,2)*pow(-1 + pow(mu,2),2)*(3*(3 + 8*b*dt)*mu + (-8 + 88*b*dt - 84*dt*gp)*pow(mu,3) - 8*(19 + 74*b*dt - 62*dt*gp)*pow(mu,5) + 
        8*(59 + 108*b*dt - 157*dt*gp)*pow(mu,7) - 16*(32 + 24*b*dt - 107*dt*gp)*pow(mu,9) - 96*(-2 + 13*dt*gp)*pow(mu,11) + 384*dt*gp*pow(mu,13) - 
        sqrt(1 - pow(mu,2)) - 4*(5 + 18*b*dt - 3*dt*gp)*pow(mu,2)*sqrt(1 - pow(mu,2)) + 16*(8 + 15*b*dt + 2*dt*gp)*pow(mu,4)*sqrt(1 - pow(mu,2)) - 
        8*(33 + 36*b*dt + 77*dt*gp)*pow(mu,6)*sqrt(1 - pow(mu,2)) + 16*(14 + 8*b*dt + 127*dt*gp)*pow(mu,8)*sqrt(1 - pow(mu,2)) - 
        32*(2 + 81*dt*gp)*pow(mu,10)*sqrt(1 - pow(mu,2)) + 1152*dt*gp*pow(mu,12)*sqrt(1 - pow(mu,2))) + 
     pow(gp,2)*(9*mu + 6*(5 + 4*b*dt)*pow(mu,3) + (-171 + 112*b*dt + 60*dt*gp)*pow(mu,5) - 4*(-25 + 198*b*dt + 53*dt*gp)*pow(mu,7) + 
        8*(59 + 190*b*dt + 22*dt*gp)*pow(mu,9) - 8*(119 + 156*b*dt - 27*dt*gp)*pow(mu,11) + 16*(44 + 24*b*dt - 33*dt*gp)*pow(mu,13) + 
        32*(-6 + 13*dt*gp)*pow(mu,15) - 128*dt*gp*pow(mu,17) - sqrt(1 - pow(mu,2)) - 30*pow(mu,2)*sqrt(1 - pow(mu,2)) - 
        3*(-13 + 40*b*dt + 4*dt*gp)*pow(mu,4)*sqrt(1 - pow(mu,2)) + 4*(44 + 106*b*dt - 17*dt*gp)*pow(mu,6)*sqrt(1 - pow(mu,2)) - 
        16*(32 + 37*b*dt - 22*dt*gp)*pow(mu,8)*sqrt(1 - pow(mu,2)) + 8*(69 + 52*b*dt - 37*dt*gp)*pow(mu,10)*sqrt(1 - pow(mu,2)) - 
        16*(18 + 8*b*dt + 29*dt*gp)*pow(mu,12)*sqrt(1 - pow(mu,2)) + 32*(2 + 27*dt*gp)*pow(mu,14)*sqrt(1 - pow(mu,2)) - 
        384*dt*gp*pow(mu,16)*sqrt(1 - pow(mu,2))) - 4*gm*gp*mu*(-1 + pow(mu,2))*
      (-1 - 5*(-2 + 3*dt*gp)*pow(mu,2) + (15 + 32*dt*gp)*pow(mu,4) + 2*(-86 + 31*dt*gp)*pow(mu,6) + (356 - 290*dt*gp)*pow(mu,8) + 
        4*(-76 + 107*dt*gp)*pow(mu,10) + (96 - 312*dt*gp)*pow(mu,12) + 96*dt*gp*pow(mu,14) + 3*(1 + dt*gp)*mu*sqrt(1 - pow(mu,2)) + 
        (-47 + 20*dt*gp)*pow(mu,3)*sqrt(1 - pow(mu,2)) + 2*(76 - 37*dt*gp)*pow(mu,5)*sqrt(1 - pow(mu,2)) - 
        2*(102 + 29*dt*gp)*pow(mu,7)*sqrt(1 - pow(mu,2)) + 4*(32 + 119*dt*gp)*pow(mu,9)*sqrt(1 - pow(mu,2)) - 
        8*(4 + 81*dt*gp)*pow(mu,11)*sqrt(1 - pow(mu,2)) + 288*dt*gp*pow(mu,13)*sqrt(1 - pow(mu,2)) + 
        2*b*dt*(-3 + 4*pow(mu,2))*(1 + 5*pow(mu,2) - 30*pow(mu,4) + 48*pow(mu,6) - 24*pow(mu,8) - 5*mu*sqrt(1 - pow(mu,2)) + 
           14*pow(mu,3)*sqrt(1 - pow(mu,2)) - 16*pow(mu,5)*sqrt(1 - pow(mu,2)) + 8*pow(mu,7)*sqrt(1 - pow(mu,2)))))/
   (2.*sqrt(2)*(-1 + pow(mu,2))*(gm*(-1 + pow(mu,2))*(1 + 8*pow(mu,2) - 24*pow(mu,4) + 16*pow(mu,6) - 4*mu*sqrt(1 - pow(mu,2)) + 
          8*pow(mu,5)*sqrt(1 - pow(mu,2))) + gp*(-1 - 11*pow(mu,2) + 36*pow(mu,4) - 40*pow(mu,6) + 16*pow(mu,8) + 6*mu*sqrt(1 - pow(mu,2)) - 
          4*pow(mu,3)*sqrt(1 - pow(mu,2)) - 8*pow(mu,5)*sqrt(1 - pow(mu,2)) + 8*pow(mu,7)*sqrt(1 - pow(mu,2))))))*ground_state;
}

// Rewrite the main using the correct ensemble - only Phi's modified
int main () {
  srand(time(NULL));
  ofstream out;
  out.open("driven.txt");
  out << (int)(tmax/dt) << endl << tmax << endl << dt << endl;

  Vector2cd initialState;
  initialState << .7, .8;
  initialState.normalize();

  int Nsteps = (int)((tmax)/dt) + 1;
  MatrixXd observables(Ncopies, Nsteps), coherences(Ncopies, Nsteps), errors(Ncopies, Nsteps);

  cout << "Using " << Ncopies << " copies, each with " << Nensemble << " states. dt = " << dt << endl;
  for (int Ncopy = 0; Ncopy < Ncopies; ++Ncopy) {
    cout << "Running copy " << Ncopy+1 << "/" << Ncopies << "...\n";

    Vector2cd psi = initialState.normalized(), psiP, psiPlus, psiMin;
    Matrix2cd exact = projector(psi), rho;
    int Nstep = 0;
    Npsi = Nensemble; Nperp = 0; Nplus = 0; Nmin = 0;
    double mu = real(psi(1));
    psiP << -mu, sqrt(1.-mu*mu);
    psiP.normalize();
    psiPlus = ((psi+psiP)/sqrt(2)).normalized();
    psiMin = ((psi-psiP)/sqrt(2)).normalized();

    for (double t = 0.; t < tmax; t += dt) {
      // Old populations
      Npsi_old = Npsi; Nperp_old = Nperp; Nplus_old = Nplus; Nmin_old = Nmin;
      // Fractions of states
      double fpsi = (double)Npsi/(double)Nensemble, fperp = (double)Nperp/(double)Nensemble, fplus = (double)Nplus/(double)Nensemble, fmin = (double)Nmin/(double)Nensemble;

      rho = fpsi*projector(psi) + fperp*projector(psiP) + fplus*projector(psiPlus) + fmin*projector(psiMin);
      if (Ncopy == 0) {
        out << observable(exact) << " " << observable(projector(psi)) << endl;
        out << coherence(exact) << " " << coherence(projector(psi)) << endl;
        out << fpsi << " " << fperp << " " << fplus << " " << fmin << endl;
        out << gamma_p(t) << " " << gamma_m(t) << " " << beta(t) << endl;
      }
      observables(Ncopy, Nstep) = observable(rho);
      coherences(Ncopy, Nstep) = coherence(rho);
      errors(Ncopy, Nstep) = TD(rho, exact);

      // Updating exact
      exact += RK4(exact,t);

      double lpsi_plus, lpsi_min, lpsiP_plus, lpsiP_min, lmin_psi, lmin_psiP, lplus_psi, lplus_psiP;
      double psi_plus = 0., psi_min = 0., psiP_plus = 0., psiP_min = 0., min_psi = 0., min_psiP = 0., plus_psi = 0., plus_psiP = 0.;

      // From psi: getting rates
      Vector2cd phi = Phi(psi,t);
      Matrix2cd R = J(projector(psi),t) + .5*(phi*psi.adjoint() + psi*phi.adjoint());
      ComplexEigenSolver<Matrix2cd> eigs;
      eigs.compute(R);
      Vector2cd eigval = eigs.eigenvalues();
      Matrix2cd eigvec = eigs.eigenvectors();
      if ( (eigvec.col(0) - psiMin).norm() < .2 ) {lpsi_min = real(eigval(0)); lpsi_plus = real(eigval(1));}
      else {lpsi_min = real(eigval(1)); lpsi_plus = real(eigval(0));}
      // Checking positivity of the rates
      if (lpsi_min >= 0.) psi_min += dt*lpsi_min;
      else {min_psi -= rev_jump(lpsi_min,Nmin,Npsi)*dt; plus_psi -= rev_jump(lpsi_min,Nplus,Npsi)*dt;}
      if (lpsi_plus >= 0.) psi_plus += dt*lpsi_plus;
      else {min_psi -= rev_jump(lpsi_plus,Nmin,Npsi)*dt; plus_psi -= rev_jump(lpsi_plus,Nplus,Npsi)*dt;}

      // From psiP: getting rates
      Vector2cd phiP = PhiPerp(psiP,t);
      R = J(projector(psiP),t) + .5*(phiP*psiP.adjoint() + psiP*phiP.adjoint());
      eigs.compute(R);
      eigval = eigs.eigenvalues();
      eigvec = eigs.eigenvectors();
      if ( (eigvec.col(0) - psiMin).norm() < .2 ) {lpsiP_min = real(eigval(0)); lpsiP_plus = real(eigval(1));}
      else {lpsiP_min = real(eigval(1)); lpsiP_plus = real(eigval(0));}
      // Checking positivity of the rates
      if (lpsiP_min >= 0.) psiP_min += dt*lpsiP_min;
      else {min_psiP -= rev_jump(lpsiP_min,Nmin,Nperp)*dt; plus_psiP -= rev_jump(lpsi_min,Nplus,Nperp)*dt;}
      if (lpsiP_plus >= 0.) psiP_plus += dt*lpsiP_plus;
      else {min_psiP -= rev_jump(lpsiP_plus,Nmin,Nperp)*dt; plus_psiP -= rev_jump(lpsi_plus,Nplus,Nperp)*dt;}

      // From psiPlus: getting rates
      Vector2cd phiPlus = PhiPlus(psiPlus,t);
      R = J(projector(psiPlus),t) + .5*(phiPlus*psiPlus.adjoint() + psiPlus*phiPlus.adjoint());
      eigs.compute(R);
      eigval = eigs.eigenvalues();
      eigvec = eigs.eigenvectors();
      if ( (eigvec.col(0) - psi).norm() < .2 ) {lplus_psi = real(eigval(0)); lplus_psiP = real(eigval(1));}
      else {lplus_psi = real(eigval(1)); lplus_psiP = real(eigval(0));}
      // Checking positivity of the rates
      if (lplus_psi >= 0.) plus_psi += dt*lplus_psi;
      else {psi_plus -= rev_jump(lplus_psi,Npsi,Nplus)*dt; psiP_plus -= rev_jump(lplus_psi,Nperp,Nplus)*dt;}
      if (lplus_psiP >= 0.) plus_psiP += dt*lplus_psiP;
      else {psi_plus -= rev_jump(lplus_psiP,Npsi,Nplus)*dt; psiP_plus -= rev_jump(lplus_psiP,Nperp,Nplus)*dt;}

      // From psiMin: getting rates
      Vector2cd phiMin = PhiPlus(psiMin,t);
      R = J(projector(psiMin),t) + .5*(phiMin*psiMin.adjoint() + psiMin*phiMin.adjoint());
      eigs.compute(R);
      eigval = eigs.eigenvalues();
      eigvec = eigs.eigenvectors();
      if ( (eigvec.col(0) - psi).norm() < .2 ) {lmin_psi = real(eigval(0)); lmin_psiP = real(eigval(1));}
      else {lmin_psi = real(eigval(1)); lmin_psiP = real(eigval(0));}
      // Checking positivity of the rates
      if (lmin_psi >= 0.) min_psi += dt*lmin_psi;
      else {psi_min -= rev_jump(lmin_psi,Npsi,Nmin)*dt; psiP_min -= rev_jump(lmin_psi,Nperp,Nmin)*dt;}
      if (lmin_psiP >= 0.) min_psiP += dt*lmin_psiP;
      else {psi_min -= rev_jump(lmin_psiP,Npsi,Nmin)*dt; psiP_min -= rev_jump(lmin_psiP,Nperp,Nmin)*dt;}

      // Printing rates and probabilities
      out << lpsi_min << " " << lpsi_plus << " " << lpsiP_min << " " << lpsiP_plus << " " << lplus_psi << " " << lplus_psiP << " " << lmin_psi << " " << lmin_psiP << endl;
      out << psi_min << " " << psi_plus << " " << psiP_min << " " << psiP_plus << " " << plus_psi << " " << plus_psiP << " " << min_psi << " " << min_psiP << endl;

      // Jumping
      for (int i = 0; i < Npsi_old; ++i)
        jump(psi_plus, psi_min, Npsi, Nplus, Nmin);
      for (int i = 0; i < Nperp_old; ++i)
        jump(psiP_plus, psiP_min, Nperp, Nplus, Nmin);
      for (int i = 0; i < Nplus_old; ++i)
        jump(plus_psi, plus_psiP, Nplus, Npsi, Nperp);
      for (int i = 0; i < Nmin_old; ++i)
        jump(min_psi, min_psiP, Nmin, Npsi, Nperp);

      // This psi is wrong (no contribution by Phi) + all the others must also be done
      // Updating psi
      Matrix2cd K = H(t) - .5*I*Gamma(t);
      psi -= I*dt*K*psi + .5*dt*phi;
      psi.normalize();
      psiP -= I*dt*K*psiP + .5*dt*phiP;
      psiP.normalize();
      psiPlus -= I*dt*K*psiPlus + .5*dt*phiPlus;
      psiPlus.normalize();
      psiMin -= I*dt*K*psiMin + .5*dt*phiMin;
      psiMin.normalize();

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