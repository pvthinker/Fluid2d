TXT = Makefile README LICENSE

PYTHON = python
F2PY   = f2py



%.so : %.f90 ; $(F2PY) -m $* -c $< --opt='-O3 -fPIC'


all:
	export F2PY=$(F2PY) ; cd core && $(MAKE)

clean:
	rm -f core/*.so core/gmg/*.so


zip:	
	cd .. ; zip -r fluid2d/fluid2d_`date '+%d_%m_%Y'`.zip fluid2d/Makefile fluid2d/README fluid2d/activate.sh fluid2d/LICENSE fluid2d/INSTALL fluid2d/requirements.txt fluid2d/setup.py fluid2d/core fluid2d/experiments/ fluid2d/docs  -x \*.so \*~ \*.pyc \*.nc \.* \*checkpoint* \*pycache* \*.rst \*.png


