.. getObsAtmo documentation master file, created by
   sphinx-quickstart on Wed Nov  1 10:40:52 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to getObsAtmo's documentation!
======================================

Version: |release|

This is the documentation for the ``getObsAtmo``, 
the emulator for air transmission in visible wavelength range. This emulator aims  at providing
an atmospheric model that can be used for atmospheric monitoring and calibration at astronomic
observation sites. It is based on the interpolation of transmission grids over
parameters generated from the ``libRadtran`` a library 
for a detailed simulation radiative transfer through the earth atmosphere. 

If this emulator has been setup for the Rubin-LSST observatory, it may be used 
for any other observation site.

This package can be downloaded from `its github repository <https://github.com/LSSTDESC/getObsAtmo/>`_ .

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   

   introduction

   quickstart

   datagrid

   rebuild

   notebooks
   
   apidocs



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
