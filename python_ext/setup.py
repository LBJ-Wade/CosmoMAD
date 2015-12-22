from distutils.core import setup, Extension

module1=Extension('py_cosmo_mad',
                  libraries=['cosmomad','gsl','gslcblas'],
                  library_dirs=['/home/damonge/lib'],
                  include_dirs=['/home/damonge/include'],
                  sources=['py_cosmo_mad.c']);

with open('README.txt') as file:
    long_description = file.read()
setup(name='PyCosmoMAD',
      version='0.0',
      description='CosmoMad\'s python extension',
      author='David Alonso',
      author_email='david.alonso@astro.ox.ac.uk',
      url='http://members.ift.uam-csic.es/dmonge/Software.html',
      license='GPL v3.0',
      long_description=long_description,
      ext_modules=[module1])
