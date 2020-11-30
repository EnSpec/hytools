.. hytools documentation master file, created by
   sphinx-quickstart on Tue Nov 24 15:15:04 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HyTools documentation
===================================
HyTools is a python library for working with imaging spectroscopy
data, with a focus on terrestrial scenes. At it's core it consists of
a series of functions for reading and writing ENVI binary images in
addition to reading NEON AOP HDF files. Built on top of
these functions are a series of higher level processing tools for data
analysis which include spectral resampling, topographic correction and
BRDF correction. Other features are currently under development and
include mask generation and MNF transformation.

We have also created a series of command line tools which string
together the processing functions and provide a more streamlined
workflow for processing images.


Contents
--------

.. toctree::
   :maxdepth: 2

   installation
   basics
   algorithms
   command_line
