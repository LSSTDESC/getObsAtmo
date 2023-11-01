Quickstart
==========

Goal
----

This package provide a fast emulator of atmospheric transmissions. These transmissions depend on the observation site,
the wavelength, the airmass, the pressure and the atmospheric content in precipitable water vapor, ozone and aerosols.

Whereas for aerosols scattering, a single component is assumed using an analytic formula, the other processes
are extracted from [libradtran](http://www.libradtran.org) simulations.

Transmission pattern for 5 observation sites are proposed: ``LSST``, ``


============= ======== ========== =========== ========= =========== ========= 
**Obs site**  **LSST**  **CTIO**    **OHP**    **PDM**   **OMK**     **OSL**  
------------- -------- ---------- ----------- --------- ----------- ---------
 pressure-hPa  731.5    774.6       937.2      710.9      600.2      1013.0 
 altitude-km   2.663    2.207       0.650      2.891      4.205         0.0 
============= ======== ========== =========== ========= =========== =========

If you don't find your observation site, choose the one wiht the closest altitude 
and provides the ground true pressure.
The transmissions will be pressure corrected for your grond pressure.         

Installation
------------

      ``
      git clone git@github.com:LSSTDESC/getObsAtmo.git
      cd getObsAtmo
      python setup.py install
      ``



About libradtran
----------------

Libradtran is an atmospheric transmission full simulation package which can be downloaded
from the [web site] http://www.libradtran.org/.

.. figure:: libradtran.png

This emulator provides interpolations from atmospheric transmissions for scattering and absorption
processes photon-air which are calculated by libradtran (version 2.0.5 for this current release).  


Details
-------

The architecture of the ``getObsAtmo`` module relies on the implementetation of 3 classes
with inheritage relationship.
The base class ``ObsAtmoGrid`` implements datafile IO and multi-dimentional interpolation
of transparencies to a 

.. tabs::

   .. tab::  ``Class ObsAtmo``

      Top-leve User-Endpoint class to call the emulator. This class inherit from the implemented class ``ObsAtmoPressure``.
      This interface class should remain as stable as possible given it does not implement behaviors than can
      be improved in later versions.

   .. tab:: ``Class ObsAtmoPressure``

      This class inherit from the base class ``ObsAtmoGrid``.
      This class implement corrections due to pressure effect on atmospheric transmissions.

   .. tab:: ``Class ObsAtmoGrid``

      This is the base class which implement the interpolation of atmospheric transmissions from a series
      of dataset extracted from [libradtran](http://www.libradtran.org/).
 

Usage
-----

The access to the emulator is done as follow.
These are detailed in :doc:`apidocs`.

.. code::

   >>> from getObsAtmo.getObsAtmo import ObsAtmo
   >>> emul =  ObsAtmo()
   >>> # or
   >>> emul =  ObsAtmo('CTIO')
   >>> # or 
   >>> emul =  ObsAtmo('LSST',743.0)
   
   >>> wl = [400.,800.,900.] # define the wavelength array
   >>> am=1.2  # set the airmass
   >>> pwv =4.0  # set the precipitable water vapor in mm
   >>> oz=300. # set the ozone depth on DU
   >>> transm = emul.GetAllTransparencies(wl,am,pwv,oz)
