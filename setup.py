from setuptools import setup, find_packages
from hytools import __version__
setup(
    name='hy-tools',
    description= 'HyTools: Hyperspectral image processing library',
    version= __version__,
    license='GNU General Public License v3.0',
    url='https://github.com/EnSpec/hytools',
    author = 'Adam Chlus, Zhiwei Ye, Ting Zheng, Natalie Queally, Evan Greenberg and Philip Townsend',
    packages=find_packages(),
    install_requires=['h5py',
                      'matplotlib',
                      'numpy',
                      'pandas',
                      'ray',
                      'scikit-learn',
                      'scipy'],
    python_requires='>=3.6, !=3.9.*'
    )
