"""
Compile Fluid2d Fortran routines into Python modules
using f2py

March 2023: the compilation now by-pass 'make' and 'Makefile'

"""
import os
import glob
import subprocess

FFLAGS = "-O3 -march=native -fPIC"

srcs = ["core/fortran_advection.f90",
        "core/fortran_diag.f90",
        "core/fortran_fluxes.f90",
        "core/fortran_operators.f90",
        "core/gmg/fortran_multigrid.f90"]


pwd = os.getcwd()
null = subprocess.PIPE

print("Compile Fortran routines to Python modules")

for filename in srcs:

    src = os.path.basename(filename)
    direc = os.path.dirname(filename)
    lib = os.path.splitext(src)[0]

    localdir = f"{pwd}/{direc}"

    cmd = ["f2py", f'--opt={FFLAGS}', "-m", lib, "-c", src]

    print(" ".join(cmd))

    os.chdir(localdir)
    subprocess.call(cmd, stdout=null, stderr=null)
