from setuptools import setup
import os


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(name='TESPy',
      version='0.1.1',
      description='Thermal Engineering Systems in Python (TESPy)',
      url='http://github.com/oemof/tespy',
      author='Francesco Witte',
      author_email='francesco.witte@web.de',
      long_description=read('README.rst'),
      license='GPL-3.0',
      packages=['tespy', 'tespy.components', 'tespy.tools'],
      python_requires='>=3',
      install_requires=['CoolProp >= 6.0.0',
                        'matplotlib >= 2.0.2',
                        'numpy >= 1.13.3',
                        'pandas >= 0.19.2',
                        'scipy >= 0.19.1',
                        'tabulate >= 0.8.2'])
