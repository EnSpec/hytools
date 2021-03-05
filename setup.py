from setuptools import setup, find_packages

setup(
    name='hy-tools',
    description= 'HyTools: Hyperspectral image processing library',
    version='1.1.0',
    license='MIT',
    url='https://github.com/EnSpec/hytools',
    author = 'Adam Chlus, Zhiwei Ye, Ting Zheng, Natalie Queally and Philip Townsend',
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


