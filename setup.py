from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from numpy import get_include

radx_ext = Extension(
    "radx",
    sources=["radx.pyx"],
    libraries=['NcfMdv', 'Radx', 'Mdv', 'dsserver', 'didss', 'rapformats',
               'toolsa', 'euclid', 'tdrp', 'dataport', 'netcdf_c++', 'netcdf',
               'hdf5_hl', 'hdf5', 'z', 'udunits2', 'bz2'],
    include_dirs=[get_include(), '/usr/local/include/Radx/'],
    library_dirs=['/usr/local/lib'],
    language="c++")

setup(
    name='pyRadx',
    ext_modules=[radx_ext],
    cmdclass={'build_ext': build_ext},
)
