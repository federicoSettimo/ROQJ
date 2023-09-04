roqj.o: roqj.h roqj.cpp
	g++ roqj.cpp -c -o roqj.o -std=c++20 -O3 -ffast-math -fno-math-errno

ph_cov: roqj.o Examples/ph_cov.cpp
	g++ Examples/ph_cov.cpp roqj.o -o Examples/ph_cov.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov.x
	python3 Examples/plot.py "Phase covariant, undriven" "$$\rho_{01}(t)$$" #ph_cov.png

ph_cov_non_P: roqj.o roqj_pop.o Examples/ph_cov_non_P.cpp
	g++ Examples/ph_cov_non_P.cpp roqj.o roqj_pop.o -o Examples/ph_cov_non_P.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_non_P.x
	python3 Examples/plot.py "Phase covariant non-P divisible" "$$ tr[\rho(t)\sigma_z] $$" ph_cov_non_P.png

ph_cov_nM: roqj.o Examples/ph_cov_nM.cpp
	g++ Examples/ph_cov_nM.cpp roqj.o -o Examples/ph_cov_nM.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_nM.x
	python3 Examples/plot.py "Non-P divisible" "$$\rho_{01}(t)$$" #ph_cov_nM.png

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
	python3 Examples/plot.py "Phase covariant, eigenstates of $$\sigma_{x}$$" "$$ tr[\rho(t)\sigma_{x}]$$" #Ph_cov_x.png

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