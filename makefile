roqj.o: roqj.h roqj.cpp
	g++ roqj.cpp -c -o roqj.o -std=c++20 -O3 -ffast-math -fno-math-errno


ph_cov: roqj.o Examples/ph_cov.cpp
	g++ Examples/ph_cov.cpp roqj.o -o Examples/ph_cov.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov.x
	python3 Examples/plot.py "Phase covariant, undriven" "$$\rho_{01}(t)$$" #ph_cov.png

ph_cov_x: roqj.o Examples/ph_cov_x.cpp
	g++ Examples/ph_cov_x.cpp roqj.o -o Examples/ph_cov_x.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_x.x
	python3 Examples/plot.py "Phase covariant, eigenstates of $$\sigma_{x}$$" "$$ tr[\rho(t)\sigma_{x}]$$" Ph_cov_x.png

ph_cov_nonP: roqj.o Examples/ph_cov_nonP.cpp
	g++ Examples/ph_cov_nonP.cpp roqj.o -o Examples/ph_cov_nonP.x -std=c++20 -O3 -ffast-math -fno-math-errno
	./Examples/ph_cov_nonP.x
	python3 Examples/plot.py "Phase covariant non P divisible, undriven" "$$\rho_{01}(t)$$" Ph_cov_nonP.png

ensemble: roqj.o Examples/ensemble.cpp
	g++ Examples/ensemble.cpp roqj.o -o Examples/ensemble.x -std=c++20 -O3 -ffast-math -fno-math-errno
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

roqj_arma.o: roqj_arma.h roqj_arma.cpp
	g++ roqj_arma.cpp -c -o roqj_arma.o -O2 -std=c++20 -I /Users/federico/armadillo-11.4.3/include/ -DARMA_DONT_USE_WRAPPER

ph_cov_arma: Examples/Ph_cov_eigen.cpp roqj_arma.o
	g++ Examples/ph_cov_arma.cpp roqj_arma.o -o Examples/ph_cov_arma.x -O2 -std=c++20 -I /Users/federico/armadillo-11.4.3/include/ -DARMA_DONT_USE_WRAPPER -framework Accelerate
	./Examples/ph_cov_arma.x
	python3 Examples/plot.py "Phase covariant, undriven" "$$\rho_{01}(t)$$" Ph_cov.png

clean:
	rm -f *.txt *.o *.x Examples/*.txt Examples/*.o Examples/*.x

clean_output:
	rm -f analytic.txt average.txt error.txt params.txt trajectories.txt