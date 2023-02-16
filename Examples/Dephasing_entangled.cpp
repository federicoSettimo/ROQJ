// 2 entangled qubits, dephasing of one. Eval if the collapse of the second into |g,e> causes the same on the first one
#include "../roqj.h"
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;
using namespace Eigen;

const double gamma_z = 5.;

// Partial trace over the qubit 1 and 2
MatrixXcd tr_1(const MatrixXcd &rho) {
	MatrixXcd A {{rho(0,0) + rho(2,2), rho(0,1) + rho(2,3)}, {rho(1,0) + rho(3,2), rho(1,1) + rho(3,3)}};
	return A;
}

MatrixXcd tr_2(const MatrixXcd &rho) {
	MatrixXcd A {{rho(0,0) + rho(1,1), rho(0,2) + rho(1,3)}, {rho(2,0) + rho(3,1), rho(2,2) + rho(3,3)}};
  return A;
}

// Trace distance
double TD (const MatrixXcd &rho, const MatrixXcd &sigma) {
  MatrixXcd A = rho - sigma;
  return 0.5*real((A.transpose()*A).sqrt().trace());
}

// From the single qubit case we know that rho_2 = |e> or |g>.
// If also the first qubit is in the same state according to the first one, obs -> 0
double observable (const MatrixXcd &rho) {
  return TD(tr_1(rho), tr_2(rho));
  //return real(trace(sigma_z*tr_1(rho)));
  //return real(trace(sigma_z*tr_2(rho)));
}

// Free-evolution effective Hamiltonian
MatrixXcd H (double t) {
  return MatrixXcd::Zero(4,4);
}

// J_t(rho)
MatrixXcd J (const MatrixXcd &rho, double t) {
  return gamma_z*kroneckerProduct(id, sigma_z)*rho*kroneckerProduct(id, sigma_z);
}

// Gamma(t)
MatrixXcd Gamma (double t) {
  return kroneckerProduct(id, gamma_z*id);
}

// C(t)
MatrixXcd C (const MatrixXcd &rho, double t) {
  return kroneckerProduct(id, gamma_z*id);
}


int main() {
  double tmin = 0., tmax = 10, dt = 0.01;
  int N_ensemble = 100, Ncopies = 1, dimH = 4, Ntraj = 7;
  bool printTraj = true;

  roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, dimH, printTraj, Ntraj);

  jump.set_N_traj_print(Ntraj);

  VectorXcd initialState(4);
  initialState << 0.,1./sqrt(2),1./sqrt(2.),0.;
  //initialState << 1./sqrt(2),0.,0.,1./sqrt(2.);
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}