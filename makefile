roqj.o: roqj.h roqj.cpp
	g++ roqj.cpp -c -o roqj.o -std=c++20 -O3 -ffast-math -fno-math-errno

ph_cov: roqj.o Examples/ph_cov.cpp
	g++ Examples/ph_cov.cpp roqj.o -o Examples/ph_cov.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov.x
	python3 Examples/plot.py "Eternally non-Markovian, driven" "$$ tr[\rho\sigma_z] $$" #ph_cov.png

ph_cov_non_P: roqj.o roqj_pop.o Examples/ph_cov_non_P.cpp
	g++ Examples/ph_cov_non_P.cpp roqj.o roqj_pop.o -o Examples/ph_cov_non_P.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_non_P.x
	python3 Examples/plot.py "Phase covariant non-P divisible" "$$ tr[\rho(t)\sigma_z] $$" ph_cov_non_P.png

ph_cov_nM: roqj.o Examples/ph_cov_nM.cpp
	g++ Examples/ph_cov_nM.cpp roqj.o -o Examples/ph_cov_nM.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_nM.x
	python3.10 Examples/plot.py "Non-P divisible" "$$\rho_{01}(t)$$" #ph_cov_nM.png

ph_cov_nM_fixed_jump: roqj.o Examples/ph_cov_nM_fixed_jump.cpp
	g++ Examples/ph_cov_nM_fixed_jump.cpp roqj.o -o Examples/ph_cov_nM_fixed_jump.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_nM_fixed_jump.x
	python3.10 Examples/plot.py "Non-P divisible" "$$\rho_{01}(t)$$" #ph_cov_nM.png

gamma_p_nM: roqj.o Examples/gamma_p_nM.cpp
	g++ Examples/gamma_p_nM.cpp roqj.o -o Examples/gamma_p_nM.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/gamma_p_nM.x
	python3.10 Examples/plot_gamma_p.py "$$\gamma_{+}<0$$" "$$ tr[\rho \sigma_z]$$"

det_evol_gamma_p: Examples/det_evol_gamma_p.cpp
	g++ Examples/det_evol_gamma_p.cpp -o Examples/det_evol_gamma_p.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/det_evol_gamma_p.x
	python3.10 Examples/det_evol_gamma_p.py

det_evol_gamma_pm: Examples/det_evol_gamma_pm.cpp
	g++ Examples/det_evol_gamma_pm.cpp -o Examples/det_evol_gamma_pm.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/det_evol_gamma_pm.x
	python3.10 Examples/det_evol_gamma_pm.py

gamma_pm_nM: roqj.o Examples/gamma_pm_nM.cpp
	g++ Examples/gamma_pm_nM.cpp roqj.o -o Examples/gamma_pm_nM.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/gamma_pm_nM.x
	python3.10 Examples/plot_gamma_pm.py "$$\gamma_{\pm}<0$$" "$$ tr[\rho \sigma_z]$$"

gamma_pm_2states: roqj.o Examples/gamma_pm_2states.cpp
	g++ Examples/gamma_pm_2states.cpp roqj.o -o Examples/gamma_pm_2states.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/gamma_pm_2states.x
	python3.10 Examples/plot_gamma_pm.py "$$\gamma_{\pm}<0$$" "$$ tr[\rho \sigma_x]$$"

plot_nM_gp: Examples/plot_nM_gp.cpp
	g++ Examples/plot_nM_gp.cpp -o Examples/plot_nM_gp.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/plot_nM_gp.x
	python3.10 Examples/plot_nM_gp.py

target_state_unravelling: Examples/target_state_unravelling.cpp
	g++ Examples/target_state_unravelling.cpp -o Examples/target_state_unravelling.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/target_state_unravelling.x
	python3.10 Examples/target_state_unravelling.py

ph_cov_nM_state: Examples/ph_cov_nM_state.cpp roqj_state.o
	g++ Examples/ph_cov_nM_state.cpp roqj_state.o -o Examples/ph_cov_nM_state.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_nM_state.x
	python3.10 Examples/plot.py

ph_cov_nM_gpm_state: Examples/ph_cov_nM_gpm_state.cpp roqj_state.o
	g++ Examples/ph_cov_nM_gpm_state.cpp roqj_state.o -o Examples/ph_cov_nM_gpm_state.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_nM_gpm_state.x
	python3.10 Examples/plot.py

ph_cov_1bit: Examples/ph_cov_1bit.cpp roqj.o
	g++ Examples/ph_cov_1bit.cpp roqj.o -o Examples/ph_cov_1bit.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_1bit.x
	python3.10 Examples/plot.py

deph_d_dim: roqj_multiple_obs.o Examples/deph_d_dim.cpp
	g++ Examples/deph_d_dim.cpp roqj_multiple_obs.o -o Examples/deph_d_dim.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/deph_d_dim.x
	python3.10 Examples/plot_multiple_obs.py

det_evol_state: Examples/det_evol_state.cpp
	g++ Examples/det_evol_state.cpp -o Examples/det_evol_state.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/det_evol_state.x
	python3.10 Examples/det_evol_state.py

det_evol_eigenvalues: Examples/det_evol_eigenvalues.cpp
	g++ Examples/det_evol_eigenvalues.cpp -o Examples/det_evol_eigenvalues.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/det_evol_eigenvalues.x
	python3.10 Examples/det_evol_eigenvalues.py

2qubits_separable: Examples/2qubits_separable.cpp roqj_state.o
	g++ Examples/2qubits_separable.cpp roqj_state.o -o Examples/2qubits_separable.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/2qubits_separable.x
	python3.10 Examples/plot.py

2qubits_separable1: Examples/2qubits_separable1.cpp roqj_state.o
	g++ Examples/2qubits_separable1.cpp roqj_state.o -o Examples/2qubits_separable1.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/2qubits_separable1.x
	python3.10 Examples/plot.py	

driven_gamma_minus: Examples/driven_gamma_minus.cpp roqj_state.o
	g++ Examples/driven_gamma_minus.cpp roqj_state.o -o Examples/driven_gamma_minus.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_gamma_minus.x
	python3.10 Examples/plot.py	

driven_gamma_minus_ens: Examples/driven_gamma_minus_ens.cpp roqj_state.o
	g++ Examples/driven_gamma_minus_ens.cpp roqj_state.o -o Examples/driven_gamma_minus_ens.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_gamma_minus_ens.x
	python3.10 Examples/driven_gamma_minus_ens.py	

driven_gamma_minus_ens_new: Examples/driven_gamma_minus_ens_new.cpp roqj_state.o
	g++ Examples/driven_gamma_minus_ens_new.cpp roqj_state.o -o Examples/driven_gamma_minus_ens_new.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_gamma_minus_ens_new.x
	python3.10 Examples/driven_gamma_minus_ens.py	

driven: Examples/driven.cpp
	g++ Examples/driven.cpp -o Examples/driven.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven.x
	python3.10 Examples/driven.py	

driven_2: Examples/driven_2.cpp
	g++ Examples/driven_2.cpp -o Examples/driven_2.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_2.x
	python3.10 Examples/driven_2.py	

driven_2ens: Examples/driven_2ens.cpp
	g++ Examples/driven_2ens.cpp -o Examples/driven_2ens.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_2ens.x
	python3.10 Examples/driven_2ens.py	

driven_2ens_gpm_equal: Examples/driven_2ens_gpm_equal.cpp
	g++ Examples/driven_2ens_gpm_equal.cpp -o Examples/driven_2ens_gpm_equal.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_2ens_gpm_equal.x
	python3.10 Examples/driven_2ens_gpm_equal.py	

driven_error: Examples/driven_error.cpp
	g++ Examples/driven_error.cpp -o Examples/driven_error.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_error.x
	python3.10 Examples/driven_error.py

enm_pm: roqj_state.o Examples/enm_pm.cpp
	g++ Examples/enm_pm.cpp roqj_state.o -o Examples/enm_pm.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/enm_pm.x
	python3.10 Examples/plot.py "Eternally non-Markovian, jumps to $$ |\pm> $$ " "$$ tr[\rho\sigma_z] $$" #ph_cov.png

driven_adaptive_dt: Examples/driven_adaptive_dt.cpp
	g++ Examples/driven_adaptive_dt.cpp -o Examples/driven_adaptive_dt.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_adaptive_dt.x
	python3.10 Examples/driven_adaptive_dt.py

driven_adaptive_4d: Examples/driven_adaptive_4d.cpp
	g++ Examples/driven_adaptive_4d.cpp -o Examples/driven_adaptive_4d.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_adaptive_4d.x
	python3.10 Examples/driven_adaptive_4d.py

driven_long_double: Examples/driven_long_double.cpp
	g++ Examples/driven_long_double.cpp -o Examples/driven_long_double.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_long_double.x
	python3.10 Examples/driven_2.py	

driven_5d_gen: Examples/driven_5d_gen.cpp
	g++ Examples/driven_5d_gen.cpp -o Examples/driven_5d_gen.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_5d_gen.x
	python3.10 Examples/driven_5d_gen.py

driven_5d_rectangle: Examples/driven_5d_rectangle.cpp
	g++ Examples/driven_5d_rectangle.cpp -o Examples/driven_5d_rectangle.x -std=c++20 -O3
	./Examples/driven_5d_rectangle.x
	python3.10 Examples/driven_5d_rectangle.py

driven_roqj_state: Examples/driven_roqj_state.cpp roqj_state.o
	g++ Examples/driven_roqj_state.cpp roqj_state.o -o Examples/driven_roqj_state.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_roqj_state.x
	python3.10 Examples/plot.py "Strong driving: $$ \beta = 10 \gamma_\pm $$" "$$ tr[\rho\sigma_z] $$"

driven_roqj_nM: Examples/driven_roqj_nM.cpp roqj_state.o
	g++ Examples/driven_roqj_nM.cpp roqj_state.o -o Examples/driven_roqj_nM.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_roqj_nM.x
	python3.10 Examples/plot.py "Driven, $$ \gamma_-<0 $$ from $$ t = \pi/2$$ " "$$ tr[\rho\sigma_z] $$"

driven_nM_no_det_evol: Examples/driven_nM_no_det_evol.cpp roqj_state.o
	g++ Examples/driven_nM_no_det_evol.cpp roqj_state.o -o Examples/driven_nM_no_det_evol.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_nM_no_det_evol.x
	python3.10 Examples/plot.py "Driven, $$ \gamma_-<0 $$ from $$ t = \pi/2$$ " "$$ tr[\rho\sigma_z] $$"

enm_driven_5d_ens: Examples/enm_driven_5d_ens.cpp roqj_state.o
	g++ Examples/enm_driven_5d_ens.cpp roqj_state.o -o Examples/enm_driven_5d_ens.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/enm_driven_5d_ens.x
	python3.10 Examples/plot.py "EnM driven, $$ 2\beta = 1- \tanh t $$" "$$ tr[\rho\sigma_z] $$"

3d_ph_cov_fixed_basis: Examples/3d_ph_cov_fixed_basis.cpp
	g++ Examples/3d_ph_cov_fixed_basis.cpp -o Examples/3d_ph_cov_fixed_basis.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/3d_ph_cov_fixed_basis.x
	python3.10 Examples/3d_ph_cov_fixed_basis.py

3d_ph_cov_roqj: Examples/3d_ph_cov_roqj.cpp roqj_state.o
	g++ Examples/3d_ph_cov_roqj.cpp roqj_state.o -o Examples/3d_ph_cov_roqj.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/3d_ph_cov_roqj.x
	python3.10 Examples/plot.py	"3d phase covariant, fixed basis, $$\gamma_{0,1} = -\tanh(t)$$" "$$\rho_{0,1}$$" #ph_cov.png

plot_eigs_det_5d: Examples/plot_eigs_det_5d.cpp
	g++ Examples/plot_eigs_det_5d.cpp -o Examples/plot_eigs_det_5d.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/plot_eigs_det_5d.x
	python3.10 Examples/plot_eigs_det_5d.py

driven_4d: Examples/driven_4d.cpp
	g++ Examples/driven_4d.cpp -o Examples/driven_4d.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_4d.x
	python3.10 Examples/driven_4d.py	

driven_gamma_minus_no_det_evol: Examples/driven_gamma_minus_no_det_evol.cpp roqj_state.o
	g++ Examples/driven_gamma_minus_no_det_evol.cpp roqj_state.o -o Examples/driven_gamma_minus_no_det_evol.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/driven_gamma_minus_no_det_evol.x
	python3.10 Examples/plot.py	

PD: Examples/PD.cpp
	g++ Examples/PD.cpp -o Examples/PD.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/PD.x
	python3.10 Examples/PD.py

ph_cov_orthogonal: roqj.o Examples/ph_cov_orthogonal.cpp
	g++ Examples/ph_cov_orthogonal.cpp roqj.o -o Examples/ph_cov_orthogonal.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_orthogonal.x
	python3 Examples/plot.py "Phase covariant, orthogonal jumps" "$$\rho_{01}(t)$$" #ph_cov.png

kappa_max_nonP: Examples/kappa_max_nonP.cpp
	g++ Examples/kappa_max_nonP.cpp -o Examples/kappa_max_nonP.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/kappa_max_nonP.x
	python3 Examples/kappa_max_nonP.py

eigs_non_P_ph_cov: Examples/eigs_non_P_ph_cov.cpp
	g++ Examples/eigs_non_P_ph_cov.cpp -o Examples/eigs_non_P_ph_cov.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/eigs_non_P_ph_cov.x
	python3 Examples/eigs_non_P_ph_cov.py

det_nM: Examples/det_nM.cpp
	g++ Examples/det_nM.cpp -o Examples/det_nM.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/det_nM.x
	python3 Examples/plot_det_nM.py

det_evol: Examples/det_evol.cpp
	g++ Examples/det_evol.cpp -o Examples/det_evol.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/det_evol.x
	python3 Examples/det_evol.py

3d_ph_cov: roqj.o Examples/3d_ph_cov.cpp
	g++ Examples/3d_ph_cov.cpp roqj.o -o Examples/3d_ph_cov.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/3d_ph_cov.x
	python3 Examples/plot.py "Phase covariant 3d" "$$ <1|\rho|1> - <2|\rho|2> $$" #3d_ph_cov.png

norm: Examples/norm.cpp
	g++ Examples/norm.cpp -o Examples/norm.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/norm.x
	python3 Examples/plot_norm.py

ph_cov_2_qubits: roqj.o Examples/ph_cov_2_qubits.cpp
	g++ Examples/ph_cov_2_qubits.cpp roqj.o -o Examples/ph_cov_2_qubits.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_2_qubits.x
	python3 Examples/plot.py "Phase covariant 2 qubits" "$$ H[tr_2\rho] $$" #ph_cov_2_qubits.png

ph_cov_x: roqj.o Examples/ph_cov_x.cpp
	g++ Examples/ph_cov_x.cpp roqj.o -o Examples/ph_cov_x.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_x.x
	python3.10 Examples/plot.py "Phase covariant, eigenstates of $$\sigma_{x}$$" "$$ tr[\rho(t)\sigma_{x}]$$" #Ph_cov_x.png

Lambda: roqj_state.o Examples/Lambda.cpp
	g++ Examples/Lambda.cpp roqj_state.o -o Examples/Lambda.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/Lambda.x
	python3 Examples/plot.py "Lambda model, CP divisible, $$ \Gamma_{31}\ne\Gamma_{32} $$ " "$$ <dark|\rho|dark> $$"

Lambda_basis_H: roqj_state.o Examples/Lambda_basis_H.cpp
	g++ Examples/Lambda_basis_H.cpp roqj_state.o -o Examples/Lambda_basis_H.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/Lambda_basis_H.x
	python3 Examples/plot.py "Lambda model, CP divisible" "$$ <dark|\rho|dark> $$"

Lambda+V_double_channel: Examples/Lambda+V_double_channel.cpp
	g++ Examples/Lambda+V_double_channel.cpp -o Examples/Lambda+V_double_channel.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/Lambda+V_double_channel.x
	python3 Examples/plot.py "Lambda(+V) model, CP divisible, jumps to |dark>, $$ \Gamma_{10}=\Gamma_{20} = .3 $$, $$ \Gamma_{01}=\Gamma_{01} = 1 $$ " "$$ <dark|\rho|dark> $$"

Lambda_force_dark: Examples/Lambda_force_dark.cpp
	g++ Examples/Lambda_force_dark.cpp -o Examples/Lambda_force_dark.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/Lambda_force_dark.x
	python3 Examples/plot.py "Lambda model, CP divisible, jumps to |dark>, $$ \Gamma_{31}=\Gamma_{32} $$, large driving " "$$ <dark|\rho|dark> $$"

V_system: roqj_state.o Examples/V_system.cpp
	g++ Examples/V_system.cpp roqj_state.o -o Examples/V_system.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/V_system.x
	python3 Examples/plot.py "V model, CP divisible" "$$ <0|\rho|0> $$"

ensemble: roqj.o roqj_mixed.o Examples/ensemble.cpp
	g++ Examples/ensemble.cpp roqj.o roqj_mixed.o -o Examples/ensemble.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ensemble.x
	python3 Examples/plot.py "Phase covariant, using an initial ensemble" "$$\rho_{01}(t)$$" #roqj_ensemble.png

ph_cov_no_plot: roqj.o Examples/ph_cov.cpp
	g++ Examples/ph_cov.cpp roqj.o -o Examples/ph_cov.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov.x

ph_cov_x_no_plot: roqj.o Examples/ph_cov_x.cpp
	g++ Examples/ph_cov_x.cpp roqj.o -o Examples/ph_cov_x.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_x.x

gen_qubit: roqj.o roqj_gen_qubit.o Examples/gen_qubit.cpp
	g++ Examples/gen_qubit.cpp roqj.o roqj_gen_qubit.o -o Examples/gen_qubit.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/gen_qubit.x
	python3 Examples/plot.py "Generic qubit dynamics" "$$\rho_{01}(t)$$" #gen_qubit.png

dephasing: roqj.o Examples/Dephasing.cpp
	g++ Examples/Dephasing.cpp roqj.o -o Examples/Dephasing.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/Dephasing.x
	python3 Examples/plot.py "Dephasing" "$$ tr[\rho(t) \sigma_z] $$" #Dephasing.png

dephasing_entangled: roqj.o Examples/Dephasing_entangled.cpp
	g++ Examples/Dephasing_entangled.cpp roqj.o -o Examples/Dephasing_entangled.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/Dephasing_entangled.x
	python3 Examples/plot.py "Dephasing on 1, state $$ (|01>+|10>)/\sqrt{2} $$" "$$ D(tr_1\rho(t), tr_2\rho(t)) $$" #Dephasing_entangled.png

no_jump: Examples/no_jump.cpp
	g++ Examples/no_jump.cpp -o Examples/no_jump.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/no_jump.x
	python3 Examples/plot_noJumps.py

roqj_pop.o: roqj_pop.cpp roqj_pop.h
	g++ roqj_pop.cpp -c -o roqj_pop.o -std=c++20 -O3 -ffast-math -fno-math-errno

roqj_mixed.o: roqj_mixed.h roqj_mixed.cpp
	g++ roqj_mixed.cpp -c -o roqj_mixed.o -std=c++20 -O3 -ffast-math -fno-math-errno

roqj_gen_qubit.o: roqj_gen_qubit.cpp roqj_gen_qubit.h
	g++ roqj_gen_qubit.cpp -c -o roqj_gen_qubit.o -std=c++20 -O3 -ffast-math -fno-math-errno

roqj_state.o: roqj_state.h roqj_state.cpp
	g++ roqj_state.cpp -c -o roqj_state.o -std=c++20 -O3 -ffast-math -fno-math-errno

functions.o: functions.cpp functions.h
	g++ functions.cpp -c -o functions.o -std=c++20 -O3 -ffast-math -fno-math-errno

roqj_multiple_obs.o: roqj_multiple_obs.h roqj_multiple_obs.cpp
	g++ roqj_multiple_obs.cpp -c -o roqj_multiple_obs.o -std=c++20 -O3 -ffast-math -fno-math-errno

ph_cov_pop: Examples/ph_cov_pop.cpp roqj.o roqj_pop.o
	g++ Examples/ph_cov_pop.cpp roqj.o roqj_pop.o -o Examples/ph_cov_pop.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_pop.x
	python3 Examples/plot.py "Phase covariant, using effective populations" "$$ tr[\rho(t)\sigma_{z}]$$" #Ph_cov_pop.png

ph_cov_pop_no_plot: Examples/ph_cov_pop.cpp roqj.o roqj_pop.o
	g++ Examples/ph_cov_pop.cpp roqj.o roqj_pop.o -o Examples/ph_cov_pop.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_pop.x

roqj_arma.o: roqj_arma.h roqj_arma.cpp
	g++ roqj_arma.cpp -c -o roqj_arma.o -O2 -std=c++20 -I /Users/federico/armadillo-11.4.3/include/ -DARMA_DONT_USE_WRAPPER

ph_cov_arma: Examples/Ph_cov_eigen.cpp roqj_arma.o
	g++ Examples/ph_cov_arma.cpp roqj_arma.o -o Examples/ph_cov_arma.x -O2 -std=c++20 -I /Users/federico/armadillo-11.4.3/include/ -DARMA_DONT_USE_WRAPPER -framework Accelerate
	./Examples/ph_cov_arma.x
	python3 Examples/plot.py "Phase covariant, undriven" "$$\rho_{01}(t)$$" Ph_cov.png

clean:
	rm -f *.txt *.o *.x Examples/*.txt Examples/*.o Examples/*.x

clean_output:
	rm -f *.txt