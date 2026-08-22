.. _install:

=============
 Installation
=============


Dependencies
============

   * numpy
   * h5py
   * ray
   * scipy
   * h5netcdf
   * matplotlib
   * pandas
   * scikit-learn


Installing HyTools
==================

If Git is installed
-------------------
To download to hytools simply clone the github repository

.. code-block:: shell
     
   $ git clone https://github.com/EnSpec/hytools.git
 

and run the following command inside the hytools folder to install

.. code-block:: shell
     
   $ pip install .


or 

.. code-block:: shell
   
  $ python -m pip install git+https://github.com/EnSpec/hytools.git

or 

.. code-block:: shell

  $ uv pip install git+https://github.com/EnSpec/hytools.git



If Git is NOT installed
-----------------------

or 

Install with pip run:

.. code-block:: shell

  $ pip install hy-tools

or

.. code-block:: shell
   
  $ pip install https://github.com/EnSpec/hytools/archive/refs/heads/master.zip

or download source code and run the following command inside 
the hytools folder to install with uv

.. code-block:: shell

  $ uv pip install .