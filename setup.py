from setuptools import setup, find_packages

setup(
    name='hy-tools',
    description= 'HyTools: Hyperspectral image processing library',
    version= '1.6.0',
    license='GNU General Public License v3.0',
    url='https://github.com/EnSpec/hytools',
    author = 'Adam Chlus, Zhiwei Ye, Ting Zheng, Natalie Queally, Evan Greenberg and Philip Townsend',
    packages=find_packages(),
    install_requires=['h5py',
                      'h5netcdf',
                      'matplotlib',
                      'numpy',
                      'pandas',
                      'ray',
                      'scikit-learn',
                      'scipy'],
    python_requires='>3.9'
    )
