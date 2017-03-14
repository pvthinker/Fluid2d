%.so : %.f90 ; $(F2PY) -m $* -c $< --opt='-O3 -fPIC' 


all:	fortran_diag.so fortran_operators.so fortran_advection.so
	cd gmg && $(MAKE)
