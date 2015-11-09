from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

cmdclass = { }

ext_modules = [
    Extension(
        "pirpy.photometry.psf_kernels",
         ["pirpy/photometry/psf_kernels.pyx"],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp']
        ),
        ]

cmdclass.update({ 'build_ext': build_ext })

setup(
    name='pirpy',
    cmdclass = cmdclass,
    ext_modules=ext_modules
)
