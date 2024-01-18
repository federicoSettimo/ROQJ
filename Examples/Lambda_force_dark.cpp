/*
    Lambda model - CP divisible
    Forcing jumps to the dark state (degenerate subspace of J for Gamma31 = Gamma32)
    Omega1 = Omega2 = 1
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

static complex<double> I(0,1);

Vector3cd three_ket {{1.,0.,0.}}, two_ket {{0.,1.,0.}}, one_ket {{0.,0.,1.}};
Vector3cd dark = (two_ket - one_ket)/sqrt(2.), phi_p = .5*(one_ket + two_ket + sqrt(2.)*three_ket), phi_m = .5*(one_ket + two_ket - sqrt(2.)*three_ket), plus_ket = (one_ket+two_ket)/sqrt(2.);
Matrix3cd sigma31 = one_ket*three_ket.adjoint(), sigma13 = three_ket*one_ket.adjoint(), sigma32 = two_ket*three_ket.adjoint(), sigma23 = three_ket*two_ket.adjoint();
double Gamma31 = 1., Gamma32 = Gamma31;

Matrix3cd comm (const Matrix3cd &A, const Matrix3cd &B) {return A*B-B*A;}

Matrix3cd anticomm (const Matrix3cd &A, const Matrix3cd &B) {return A*B+B*A;}

Matrix3cd projector (const Vector3cd &psi) {return psi*psi.adjoint();}

Matrix3cd H (double t) {
    return three_ket*two_ket.adjoint() + two_ket*three_ket.adjoint() + three_ket*one_ket.adjoint() + one_ket*three_ket.adjoint();
}

Matrix3cd Gamma (double t) {
    return Gamma31*sigma13*sigma31 + Gamma32*sigma23*sigma32;
    //return (Gamma31 + Gamma32)*projector(three_ket);
}

Matrix3cd J (const Matrix3cd &rho, double t) {
    return Gamma31*sigma31*rho*sigma13 + Gamma32*sigma32*rho*sigma23;
    //return Gamma31*one_ket*three_ket.adjoint()*rho*three_ket*one_ket.adjoint() + Gamma32*two_ket*three_ket.adjoint()*rho*three_ket*two_ket.adjoint();
}

Vector3cd Phi(const Vector3cd &psi, double t) {// Phi for after the jump (not to dark)
    return 0.*three_ket;
}

double observable (const Matrix3cd &rho) {
    return real((rho*projector(dark)).trace());
}

int main () {
    double tmin = 0., tmax = 5, dt = 0.001;
    int N_ensemble = 1000, Ncopies = 3, dimH = 3, Ntraj = 5;
    bool printTraj = true;

    ofstream params;
    params.open("params.txt");
    params << 1 << endl << N_ensemble << endl << tmin << endl << tmax << endl << dt << endl << true << endl << Ntraj << endl << 3;
    params.close();

    Vector3cd initialState {{1.,0.,1.}}; initialState.normalize();
    std::vector<Vector3cd> psi;
    vector<bool> jumped(N_ensemble);
    for (int i = 0; i <= N_ensemble; ++i) {
        psi.push_back(initialState);
        jumped[i] = false;
    }

    Matrix3cd rho_ex = projector(initialState);
    ofstream out_ex, traj, out_obs, out_err;
    out_ex.open("analytic.txt");
    traj.open("trajectories.txt");
    out_obs.open("average.txt");
    out_err.open("error.txt");

    for (double t = tmin; t <= tmax; t += dt) {
        // Prints and evolves the exact solution
        out_ex << observable(rho_ex) << endl;
        rho_ex = rho_ex + (-complex<double>(0,1)*comm(H(t),rho_ex) + J(rho_ex,t) - 0.5*anticomm(Gamma(t),rho_ex))*dt;
        rho_ex /= rho_ex.trace();

        // Average state
        Matrix3cd rho;

        // Cycle on the ensemble members
        for (int i = 0; i < N_ensemble; ++i) {
            if (real(psi[i](0)) < 0. || (real(psi[i](0)) < .001 && real(psi[i](1)) < 0.)) psi[i] *= -1.;
            // Prints the trajectories
            if (i < Ntraj)
                traj << observable(projector(psi[i])) << " ";

            // Updates the average
            rho += projector(psi[i])/((double)N_ensemble);
            
            
            // Draws a random number and calculates whether the evolution is deterministic or via a jump
            double z = (double)rand()/((double)RAND_MAX);

            //cout << psi[i] << endl;
            //cout << J(projector(psi[i]),t) << endl << endl;

            if (!jumped[i]) {
                complex<double> alpha3 = psi[i](0);
                if (z <= Gamma31*norm(alpha3)*dt) {
                    psi[i] = dark;
                    jumped[i] = true;
                }
                else if (z <= 2.*Gamma31*norm(alpha3)*dt) {
                    psi[i] =  (two_ket + one_ket)/sqrt(2.);
                    jumped[i] = true;
                }
                else {
                    Matrix3cd K = H(t) - 0.5*I*Gamma(t);
                    psi[i] = psi[i] - K*psi[i]*I*dt;
                }
            }
            else if (abs(dark.dot(psi[i])) < 1.) { // If in dark, nothing happens anymore
                complex<double> alpha_p = phi_p.dot(psi[i]), alpha_m = phi_m.dot(psi[i]);
                if (z <= .5*Gamma31*norm(alpha_p-alpha_m)*dt)
                    psi[i] = dark;
                else if (z <= Gamma31*norm(alpha_p-alpha_m)*dt)
                    psi[i] = (phi_m + phi_p)/sqrt(2.);
                else {
                    Matrix3cd K = H(t) - 0.5*I*Gamma(t);
                    psi[i] = psi[i] - K*psi[i]*I*dt;
                }
            }
            psi[i] = psi[i].normalized();
        }
        // Storing the observable
        out_obs << observable(rho) << endl;
        out_err << 0 << endl;

        traj << endl;
    }

    return 0;
}