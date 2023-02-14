// 2 entangled qubits, dephasing of one. Eval if the collapse of the second into |g,e> causes the same on the first one
#include "../roqj.h"

using namespace std;
using namespace arma;

const double gamma_z = 5.;

complex<double> I(0,1), one(1,0);
cx_mat sigma_z = {{one,0},{0,-one}}, Id = {{one,0},{0,one}}, sigma_x = {{0,one},{one,0}}, sigma_y = {{0,-I},{I,0}};

// Partial trace over the qubit 1 and 2
cx_mat tr_1(const cx_mat &rho) {
	cx_mat A = {{rho(0,0) + rho(2,2), rho(0,1) + rho(2,3)}, {rho(1,0) + rho(3,2), rho(1,1) + rho(3,3)}};
	return A;
}

cx_mat tr_2(const cx_mat &rho) {
	cx_mat A = {{rho(0,0) + rho(1,1), rho(0,2) + rho(1,3)}, {rho(2,0) + rho(3,1), rho(2,2) + rho(3,3)}};
  return A;
}

// Trace distance
double TD (const cx_mat &rho, const cx_mat &sigma) {
  cx_mat A = rho - sigma;
  return 0.5*real(trace( sqrtmat_sympd(A.t()*A) ));
}

// From the single qubit case we know that rho_2 = |e> or |g>.
// If also the first qubit is in the same state according to the first one, obs -> 0
double observable (const cx_mat &rho) {
  return TD(tr_1(rho), tr_2(rho));
  //return real(trace(sigma_z*tr_1(rho)));
  //return real(trace(sigma_z*tr_2(rho)));
}

// Free-evolution effective Hamiltonian
cx_mat H (double t) {
  return cx_mat(4,4,fill::zeros);
}

// J_t(rho)
cx_mat J (const cx_mat &rho, double t) {
  return 0.5*gamma_z*kron(Id, sigma_z)*rho*kron(Id, sigma_z);
}

// Gamma(t)
cx_mat Gamma (double t) {
  return kron(Id, .5*gamma_z*Id);
}

// C(t)
cx_mat C (const cx_mat &rho, double t) {
  return kron(Id, (1.-exp(-t))*gamma_z*Id);
}


int main() {
  double tmin = 0., tmax = 10, dt = 0.01;
  int N_ensemble = 100, Ncopies = 1, dimH = 4, Ntraj = 7;
  bool printTraj = true;

  roqj jump(N_ensemble, tmin, tmax, dt, Ncopies, dimH, printTraj, Ntraj);

  jump.set_N_traj_print(Ntraj);

  //cx_vec initialState = {1./sqrt(2),0.,0.,1./sqrt(2.)};
  cx_vec initialState = {0.,1./sqrt(2),1./sqrt(2.),0.};
  jump.set_initial_state(initialState);

  jump.run();

  jump.get_observable("average.txt");
  jump.get_error_observable("error.txt");

  return 0;
}