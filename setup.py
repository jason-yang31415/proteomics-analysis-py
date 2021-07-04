from setuptools import setup, find_packages

setup(
    name='proteomics-analysis',
    version='0.0.1',
    author='Jason Yang',
    author_email='jy31415@mit.edu',
    description='basic analysis for LFQ proteomics data',
    package_dir={'': './src'},
    packages=find_packages(where='./src')
)
