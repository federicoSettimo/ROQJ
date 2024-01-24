/*
    Lambda model - CP divisible
    Omega1 = Omega2 = 1
*/

#include "../roqj_state.h"

using namespace std;
using namespace Eigen;

Vector3cd three_ket {{1.,0.,0.}}, two_ket {{0.,1.,0.}}, one_ket {{0.,0.,1.}};
Vector3cd dark = (two_ket - one_ket)/sqrt(2.), phi_p = .5*(one_ket + two_ket + sqrt(2.)*three_ket), phi_m = .5*(one_ket + two_ket - sqrt(2.)*three_ket);
Matrix3cd sigma31 = one_ket*three_ket.adjoint(), sigma13 = three_ket*one_ket.adjoint(), sigma32 = two_ket*three_ket.adjoint(), sigma23 = three_ket*two_ket.adjoint();
double Gamma31 = 1., Gamma32 = 1.;

MatrixXcd H (double t) {
    return three_ket*two_ket.adjoint() + two_ket*three_ket.adjoint() + three_ket*one_ket.adjoint() + one_ket*three_ket.adjoint();
}

MatrixXcd Gamma (double t) {
    return Gamma31*sigma13*sigma31 + Gamma32*sigma23*sigma32;
    //return (Gamma31 + Gamma32)*projector(three_ket);
}

MatrixXcd J (const MatrixXcd &rho, double t) {
    return Gamma31*sigma31*rho*sigma13 + Gamma32*sigma32*rho*sigma23;
    //return Gamma31*one_ket*three_ket.adjoint()*rho*three_ket*one_ket.adjoint() + Gamma32*two_ket*three_ket.adjoint()*rho*three_ket*two_ket.adjoint();
}

VectorXcd Phi (const VectorXcd &psi, double t, bool jumped) {
    //return 0.*dark;
    complex<double> alphaP = psi(0), alphaM = psi(1), alpha0 = psi(2);
    if (!jumped)
        return (((alphaM - alphaP)*(Gamma31 - Gamma32)*(conj(alphaM) - conj(alphaP)))/(4.*sqrt(2)*conj(alpha0)))*phi_p + 
            (((alphaM - alphaP)*(Gamma31 - Gamma32)*(conj(alphaM) - conj(alphaP)))/(4.*sqrt(2)*conj(alpha0)))*phi_m;
    
    // If in dark, return 0
    if ((psi - dark).norm() < .05)
        return 0.*dark;

    // if orthogonal to dark
    return 0.*dark;
}

double observable (const MatrixXcd &rho) {
    return real((rho*projector(dark)).trace());
}

int main () {
    double tmin = 0., tmax = 5, dt = 0.001;
    int N_ensemble = 1000, Ncopies = 3, dimH = 3, Ntraj = 5;
    bool printTraj = true;

    qutrit_roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, printTraj, Ntraj);

    Vector3cd initialState {{1.,1.,1.}}; initialState.normalize();
    jump.set_initial_state(initialState);

    jump.run();

    jump.get_observable("average.txt");
    jump.get_error_observable("error.txt");
    return 0;
}