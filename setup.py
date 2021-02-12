from setuptools import setup, find_packages

setup(
    name='HyTools',
    version='1.0',
    url='https://github.com/EnSpec/hytools',
    author = 'Adam Chlus, Zhiwei Ye, Ting Zheng, Natalie Queally and Philip Townsend'
    packages=find_packages(),
    install_requires=['h5py>=3.1.0',
                      'matplotlib>=3.3.3',
                      'numpy>=1.19.4',
                      'pandas>=1.1.0',
                      'scikit-learn>=0.23.0',
                      'scipy>=1.5.3']
    )


