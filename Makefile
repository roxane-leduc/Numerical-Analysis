FC=gfortran
OPT=-g

mhouseholder: mhouseholder.o triangularization.o eigenvalues.o decompose_QR.o 
	$(FC) $(OPT) mhouseholder.o triangularization.o eigenvalues.o decompose_QR.o -o mhouseholder

mhouseholder.o: mhouseholder.f
	$(FC) $(OPT) mhouseholder.f -c

triangularization.o: triangularization.f
	$(FC) $(OPT) triangularization.f -c

eigenvalues.o: eigenvalues.f
	$(FC) $(OPT) eigenvalues.f -c

decompose_QR.o: decompose_QR.f
	$(FC) $(OPT) decompose_QR.f -c

clean:
	rm *.o mhouseholder
