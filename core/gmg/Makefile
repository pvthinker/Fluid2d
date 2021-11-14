%.so : %.f90 ; $(F2PY) --quiet -m $* -c $< --opt='-O3 -fPIC -Wno-unused-variable -Wno-maybe-uninitialized -Wno-unused-dummy-argument -Wno-tabs -Wno-cpp -Wno-unused-function'


all:	fortran_multigrid.so
