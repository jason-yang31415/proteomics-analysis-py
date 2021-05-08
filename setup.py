from setuptools import setup, find_packages

setup(
    name='protpipe-jy',
    version='0.0.1',
    author='Jason Yang',
    author_email='jy31415@mit.edu',
    description='basic proteomics pipeline to analyze LFQ data',
    package_dir={'': './src'},
    packages=find_packages(where='./src')
)
