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

Vector3cd zero_ket {{1.,0.,0.}}, two_ket {{0.,1.,0.}}, one_ket {{0.,0.,1.}};
Vector3cd dark = (two_ket - one_ket)/sqrt(2.), phi_p = .5*(one_ket + two_ket + sqrt(2.)*zero_ket), phi_m = .5*(one_ket + two_ket - sqrt(2.)*zero_ket), light = (two_ket + one_ket)/sqrt(2.);
Matrix3cd sigma01 = one_ket*zero_ket.adjoint(), sigma10 = zero_ket*one_ket.adjoint(), sigma02 = two_ket*zero_ket.adjoint(), sigma20 = zero_ket*two_ket.adjoint();
double gamma1 = 1., gamma2 = gamma1, Gamma1 = .3, Gamma2 = Gamma1;

Matrix3cd comm (const Matrix3cd &A, const Matrix3cd &B) {return A*B-B*A;}

Matrix3cd anticomm (const Matrix3cd &A, const Matrix3cd &B) {return A*B+B*A;}

Matrix3cd projector (const Vector3cd &psi) {return psi*psi.adjoint();}

MatrixXcd H (double t) {
    return zero_ket*two_ket.adjoint() + two_ket*zero_ket.adjoint() + zero_ket*one_ket.adjoint() + one_ket*zero_ket.adjoint();
}

MatrixXcd Gamma (double t) {
    return Gamma1*sigma01*sigma10 + Gamma2*sigma02*sigma20 + gamma1*sigma10*sigma01 + gamma2*sigma20*sigma02;
}

MatrixXcd J (const MatrixXcd &rho, double t) {
    return Gamma1*sigma10*rho*sigma01 + Gamma2*sigma20*rho*sigma02 + gamma1*sigma01*rho*sigma10 + gamma2*sigma02*rho*sigma20;
}

double observable (const MatrixXcd &rho) {
    return real((rho*projector(dark)).trace());
}

int main () {
    double tmin = 0., tmax = 5, dt = 0.001;
    int N_ensemble = 10000, Ncopies = 3, dimH = 3, Ntraj = 7;
    bool printTraj = true;

    ofstream params;
    params.open("params.txt");
    params << 1 << endl << N_ensemble << endl << tmin << endl << tmax << endl << dt << endl << true << endl << Ntraj << endl << 3;
    params.close();

    Vector3cd initialState {{0.,1.,1.}}; initialState.normalize();
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
                complex<double> psi0 = psi[i](0);
                if (z <= gamma2*norm(psi0)*dt) {
                    psi[i] = dark;
                    jumped[i] = true;
                }
                else if (z <= 2.*gamma2*norm(psi0)*dt) {
                    psi[i] =  (two_ket + one_ket)/sqrt(2.);
                    jumped[i] = true;
                }
                else if (z <= 2.*gamma2*norm(psi0)*dt + Gamma2*(1.-norm(psi0))*dt) {
                    psi[i] =  zero_ket;
                    jumped[i] = true;
                }
                else {
                    Matrix3cd K = H(t) - 0.5*I*Gamma(t);
                    psi[i] = psi[i] - K*psi[i]*I*dt;
                }
            }
            else if (abs(dark.dot(psi[i])) < 1.) { // Orthogonal to dark
                complex<double> psi_p = phi_p.dot(psi[i]), psi_m = phi_m.dot(psi[i]);
                if (z <= .5*gamma2*norm(psi_p-psi_m)*dt)
                    psi[i] = dark;
                else if (z <= gamma2*norm(psi_p-psi_m)*dt)
                    psi[i] = (phi_p+phi_m)/sqrt(2.);
                else if (z <= gamma2*norm(psi_p-psi_m)*dt + .5*Gamma1*norm(psi_p+psi_m)*dt)
                    psi[i] = (phi_m-phi_p)/sqrt(2.);
                else {
                    Matrix3cd K = H(t) - 0.5*I*Gamma(t);
                    psi[i] = psi[i] - K*psi[i]*I*dt;
                }
            }
            else { // Now it still can jump away from dark
                if (z <= Gamma2*dt)
                    psi[i] = (phi_m-phi_p)/sqrt(2.);
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