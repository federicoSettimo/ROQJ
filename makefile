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
	#python3 Examples/plot_nM.py "Phase covariant, undriven" "$$\rho_{01}(t)$$" #ph_cov_nM.png

ph_cov_2_qubits: roqj.o Examples/ph_cov_2_qubits.cpp
	g++ Examples/ph_cov_2_qubits.cpp roqj.o -o Examples/ph_cov_2_qubits.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_2_qubits.x
	python3 Examples/plot.py "Phase covariant 2 qubits" "$$ tr[\rho_2\sigma_x] $$" #ph_cov_2_qubits.png

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