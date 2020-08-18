import setuptools
from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

INCLUDE_DIRS=['simba3d_cpp/src',numpy.get_include()]
setup (
		name='SIMBA3D',
		version='3.0.1',
		author='Michael M. Rosenthal',
		author_email='Michael.M.Rosenthal@gmail.com',
		url='',
		keywords='chromatin',
		license='MIT',
		description='Structure Inference from Multiscale BAyes in 3 Dimensions (SIMBA3D)',
		long_description=open('README.md').read(),
		packages=setuptools.find_packages(exclude=["*.tests","*.tests.*","tests.*","tests"]),
		entry_points={'console_scripts' : [
		      #'simba3d=simba3d.simba3d_taskrun:main',
		       'simba3d=simba3d.simba3d_taskrun_SOLO:main',
		       'simba3d-result-disp=simba3d.display_results:main',
		       'simba3d-task-disp=simba3d.display_task:main',
		       'simba3d-convertion=simba3d.convertion_tool:main',
		       'simba3d-print=simba3d.print_result_summary:main',
		       'simba3d-initialize=simba3d.initialize_curve:main'
		       ]},
		install_requires=[   'numpy',
		                    'scipy'
		                    ],
		extras_require={     'ploting_tools':'matplotlib'
		                    },
		ext_modules = cythonize(Extension(
		           "simba3d.cysimba3d",                                # the extension name
		           include_dirs=INCLUDE_DIRS,
		           sources=["simba3d/cysimba3d.pyx"], #W the Cython source and
		           #library_dirs=[os.getcwd()],
		           language="c++",                        # generate and compile C++ code
		      ),annotate=True)                      
       )
