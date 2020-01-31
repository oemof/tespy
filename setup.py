from setuptools import setup
import os


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(name='TESPy',
      version='0.2.1 dev',
      description='Thermal Engineering Systems in Python (TESPy)',
      url='http://github.com/oemof/tespy',
      author='Francesco Witte',
      author_email='francesco.witte@web.de',
      long_description=read('README.rst'),
      license='MIT',
      packages=['tespy', 'tespy.components', 'tespy.data', 'tespy.networks',
                'tespy.tools'],
      include_package_data=True,
      data_files=[('tespy/data', ['tespy/data/char_lines.json',
                                  'tespy/data/char_maps.json'])],
      python_requires='>=3',
      install_requires=['CoolProp >= 6.0.0',
                        'numpy >= 1.13.3',
                        'pandas <= 1.0.0',
                        'scipy >= 0.19.1',
                        'tabulate >= 0.8.2'])
