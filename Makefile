TXT = Makefile README LICENSE

PYTHON = python
F2PY   = f2py


%.so : %.f90 ; $(F2PY) -m $* -c $< --opt='-O3 -fPIC' --quiet


all:
	export F2PY=$(F2PY) ; cd core && $(MAKE)

clean:
	rm -f core/*.so core/gmg/*.so


zip:	
	cd .. ; zip -r Fluid2d/Fluid2d_`date '+%d_%m_%Y'`.zip Fluid2d/Makefile Fluid2d/README* Fluid2d/activate.* Fluid2d/LICENSE Fluid2d/INSTALL Fluid2d/requirements.txt Fluid2d/setup.py Fluid2d/core Fluid2d/experiments/ Fluid2d/docs  -x \*.so \*~ \*.pyc \*.nc \.* \*checkpoint* \*pycache* \*.rst \*.png


