"""
"""

from setuptools import setup, find_packages
import re
from codecs import open
from os import path
import sys
try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements

here = path.abspath(path.dirname(__file__))

VERSIONFILE = "annofilt/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

if sys.version_info <= (3, 0):
    sys.stderr.write("ERROR: annofilt requires Python 3.5 " +
                     "or above...exiting.\n")
    sys.exit(1)

## parse requirements file
install_reqs = parse_requirements("requirements.txt",
                                  session=False)
# see https://stackoverflow.com/questions/62114945
try:
    requirements = [str(ir.req) for ir in install_reqs]
except:
    requirements = [str(ir.requirement) for ir in install_reqs]
setup(
    name='annofilt',
    version=verstr,
    description='Filter annotations from prokka',
    long_description="Filter out missassembled/truncated annotations",
    url='https://github.com/nickp60/annofilt',
    author='Nick Waters',
    author_email='nickp60@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='bioinformatics assembly genomics development',
    packages=['annofilt'],
    install_requires=requirements,
    include_package_data=True,
    # package_data={
    #    '': [path.join(__name__, "BugBuilder", "config_data/*")],
    # },
    entry_points={
       'console_scripts': [
           'annofilt=annofilt.annofilt:main',
           'get_complete_genomes=annofilt.get_complete_genomes:main',
           'make_annofilt_pangenome=annofilt.make_annofilt_pangenome:main',
       ],
    },
)
