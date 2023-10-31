from setuptools import setup
import os
import re


# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


version_file = os.path.join('getObsAtmo', '_version.py')
verstrline = open(version_file, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    current_version = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (version_file,))
print(f'getObsAtmo version is {current_version}')

setup(
    name='getObsAtmo',
    version=current_version,
    packages=['getObsAtmo'],
    install_requires=['numpy>1.15', 'matplotlib>3.1', 'pickle', 'scipy'],
    test_suite='nose.collector',
    tests_require=['nose'],
    package_dir={'getObsAtmo': './getObsAtmo'},
    package_data={'getObsAtmo': ['../obsatmo_data/*.pickle','../obsatmo_data/*.npy']},
    url='https://github.com/LSSTDESC/getObsAtmo',
    license='BSD',
    python_requires='>=3.7',
    author='S. Dagoret-Campagne',
    author_email='sylvie.dagoret-campagne@ijclab.in2p3.fr',
    description='',
    long_description=long_description,
    long_description_content_type='text/markdown'
)
