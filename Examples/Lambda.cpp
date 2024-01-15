/*
    Lambda model - CP divisible
    Omega1 = Omega2 = 1
*/

#include "../roqj_state.h"

using namespace std;
using namespace Eigen;

Vector3cd phi_p {{1.,0.,0.}}, phi_m {{0.,1.,0.}}, dark {{0.,0.,1.}};
Vector3cd three_ket = (phi_p - phi_m)/sqrt(2), two_ket = (phi_p + phi_m, + sqrt(2)*dark)/2., one_ket = (phi_p + phi_m - sqrt(2)*dark)/2.;
double Gamma31 = 1., Gamma32 = 1.;

MatrixXcd H (double t) {
    return sqrt(2)*projector(phi_p) - sqrt(2)*projector(phi_m);
}

MatrixXcd Gamma (double t) {
    return (Gamma31 + Gamma32)*projector(three_ket);
}

MatrixXcd J (const MatrixXcd &rho, double t) {
    return Gamma31*one_ket*three_ket.adjoint()*rho*three_ket*one_ket.adjoint() + Gamma32*two_ket*three_ket.adjoint()*rho*three_ket*two_ket.adjoint();
}

VectorXcd Phi (const VectorXcd &psi, double t, bool jumped) {
    complex<double> alphaP = psi(0), alphaM = psi(1), alpha0 = psi(2);
    if (!jumped)
        return (((alphaM - alphaP)*(Gamma31 - Gamma32)*(conj(alphaM) - conj(alphaP)))/(4.*sqrt(2)*conj(alpha0)))*phi_p + 
            (((alphaM - alphaP)*(Gamma31 - Gamma32)*(conj(alphaM) - conj(alphaP)))/(4.*sqrt(2)*conj(alpha0)))*phi_m;
    
    // If in dark, return 0
    if ((psi - dark).norm() < .05)
        return 0.*dark;

    // IF orthogonal to dark
        return 0.*dark;
}

double observable (const MatrixXcd &rho) {
    return real((rho*projector(dark)).trace());
}

int main () {
    double tmin = 0., tmax = 5, dt = 0.001;
    int N_ensemble = 10000, Ncopies = 3, dimH = 3, Ntraj = 5;
    bool printTraj = true;

    roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, dimH, printTraj, Ntraj, true, 1e-2);

    Vector3cd initialState {{1.,1.,1.}};
    jump.set_initial_state(initialState);

    jump.run();

    jump.get_observable("average.txt");
    jump.get_error_observable("error.txt");
    return 0;
}