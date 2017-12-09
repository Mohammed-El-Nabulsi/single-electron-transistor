from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='HaPPPy',
    version='0.0.1',
    description='HaPPPy package for the programming project',
    long_description=readme,
    author='Lars-Hendrik Frahm',
    author_email='lfrahm@physnet.uni-hamburg.de',
    url='https://hp.physnet.uni-hamburg.de/pfannkuche/',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)

