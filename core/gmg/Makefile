%.so : %.f90 ; $(F2PY) -m $* -c $< --opt='-O3 -fPIC' 


all:	fortran_multigrid.so
