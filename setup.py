from setuptools import setup, find_packages

setup(
    name='hytools install',
    version='0.2',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'h5py',
        'ray',
        'scipy']
    )
