from setuptools import setup, find_packages

setup(
    name='hy-tools',
    description= 'HyTools: Hyperspectral image processing library',
    version='1.0',
    license='MIT',
    url='https://github.com/EnSpec/hytools',
    download_url = 'https://github.com/EnSpec/hytools/archive/1.0.0.tar.gz',
    author = 'Adam Chlus, Zhiwei Ye, Ting Zheng, Natalie Queally and Philip Townsend',
    packages=find_packages(),
    install_requires=['h5py>=3.1.0',
                      'matplotlib>=3.3.3',
                      'numpy>=1.19.4',
                      'pandas>=1.1.0',
                      'scikit-learn>=0.23.0',
                      'scipy>=1.5.3']
    )


