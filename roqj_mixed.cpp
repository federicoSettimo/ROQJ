#include "roqj_mixed.h"
// ------------------------- ROQJ initially mixed -------------------------
roqj_mixed::roqj_mixed (int N_states, double t_i, double t_f, double dt, int N_copies, int dim_Hilbert_space, bool verbose, double threshold) {
  initialize(N_states, t_i, t_f, dt, N_copies, dim_Hilbert_space, false, 0, verbose, threshold);
}

// Setting the ensemble
void roqj_mixed::set_ensemble (const vector<pair<double, VectorXcd>> &ensemble) {
  _ensemble = ensemble;
  double sum_p = 0.;
  for (auto & elem : _ensemble) {
    if (elem.first < 0.)
      elem.first = -elem.first;
    sum_p += elem.first;
  }
  if (sum_p == 0.) {
    set_ensemble();
    return;
  }
  if (sum_p != 1.) {
    for (auto & elem : _ensemble)
      elem.first /= sum_p;
  }
}

void roqj_mixed::set_ensemble (const vector<double> &probabilities, const vector<VectorXcd> &states) {
  if (probabilities.size() != states.size()) {
    cerr << "Different number of probabilities and states -- EXIT\n";
    exit(EXIT_FAILURE);
  }
  _ensemble.clear();
  vector<pair<double, VectorXcd>> ens;
  for (int i = 0; i < probabilities.size(); ++i)
    ens.push_back(pair<double, VectorXcd>(probabilities[i], states[i]));
  set_ensemble(ens);
}

void roqj_mixed::set_ensemble () {
  VectorXcd init_state = VectorXcd::Ones(_dim_Hilbert_space).normalized();
  _ensemble.clear();
  _ensemble.push_back(pair<double, VectorXcd>(1., init_state));
}

void roqj_mixed::add_ensemble (const pair<double, VectorXcd> &state) {
  pair<double, VectorXcd> ens(state);
  if (ens.first == 0.) return;
  if (ens.first < 0.) ens.first = -ens.first;
  _ensemble.push_back(ens);

  double normalization = 0.;
  for (auto & i : _ensemble)
    normalization += i.first;
  if (normalization != 1.) {
    for (auto & i : _ensemble)
      i.first /= normalization;
  }
}

void roqj_mixed::add_ensemble (double prob, const VectorXcd &state) {
  add_ensemble(pair<double, VectorXcd>(prob, state));
}

// Runs the ROQJ for each state and takes the average state. Repeats it _N_copies time
void roqj_mixed::run () {
  if (_verbose) {
    print_info();
    print_ens();
  }

  ofstream params;
  params.open("params.txt");
  params << _N_copies << endl << _N_states << endl << _t_i << endl << _t_f << endl << _dt << endl << false << endl << 0 << endl << _dim_Hilbert_space;
  params.close();

  MatrixXd matrix_observables = MatrixXd::Zero(_num_timesteps, _N_copies);
  // Cycle on the copies
  for (int i = 0; i < _N_copies; ++i) {
    if (_verbose)
      cout << "Running copy " << i+1 << "/" << _N_copies << "...\n";
    
    // Cycle on the ensemble members'
    VectorXd obs = VectorXd::Zero(_num_timesteps);
    int n = 1;
    for (auto & elem : _ensemble) {
      set_initial_state(elem.second);
      if (_verbose)
        cout << "\tEnsemble member " << n++ << "/" << _ensemble.size() << "...\n";
      obs += elem.first * run_single_iterations(false);
    }
    if (_verbose) cout << endl;

    for (int j = 0; j < _num_timesteps; ++j)
      matrix_observables(j,i) = obs(j);
  }

  for (int i = 0; i < _num_timesteps; ++i) {
    _observable[i] = matrix_observables.row(i).mean();
    _sigma_observable[i] = sqrt((matrix_observables.row(i).array() - _observable[i]).square().sum() / (matrix_observables.row(i).size() - 1));
  }
}


VectorXd roqj_mixed::get_exact_sol (string file_out) {
  ofstream out;
  bool verbose = file_out != "";
  if (verbose)
    out.open(file_out);

  MatrixXcd rho = MatrixXcd::Zero(_dim_Hilbert_space, _dim_Hilbert_space);
  for (auto & elem : _ensemble)
    rho += elem.first*projector(elem.second);

  VectorXd obs = VectorXd(_num_timesteps);
  int n = 0;

  for (double t = _t_i; t <= _t_f; t += _dt) {
    double this_obs = observable(rho);
    if (verbose)
      out << this_obs << endl;
    obs[n] = this_obs;
    n++;
    rho = rho + (-I*comm(H(t),rho) + J(rho,t) - 0.5*anticomm(Gamma(t),rho))*_dt;
  }

  return obs;
}

void roqj_mixed::print_ens () {
  cout << "Ensemble:\n";
  cout << "\tState\t| Probability\n";
  for (auto & elem : _ensemble) {
    cout << "----------------|----------------\n";
    cout << elem.second << "\t| " << elem.first << endl;
  }
  cout << endl;
  return;
}