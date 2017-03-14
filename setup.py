#python  setup.py build_ext --inplace
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension


def gen_ext(module,fortran_options):  
    src = module.replace('.','/')+'.f90'
    ext = Extension(name=module,
                    sources=[src], # <- should be a list type
                    extra_f90_compile_args=fortran_options)
    return ext

fortran_options = ['-O3', '-fPIC']
fortran_modules = ['core.fortran_advection',
                   'core.fortran_diag',
                   'core.fortran_operators',
                   'core.gmg.fortran_multigrid']

extensions = []
for m in fortran_modules:
    extensions.append( gen_ext(m,fortran_options) )

setup(
    name='Fluid2d',
    version='1.42',
    author='Guillaume Roullet',
    author_email='roullet@univ-brest.fr',
    license='GPL3',
    packages = ['core',
                'core.gmg',
                'exp'],
    ext_modules = extensions )
