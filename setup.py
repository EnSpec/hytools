from setuptools import setup, find_packages

setup(
    name='hytools setup test',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'h5py']
    )