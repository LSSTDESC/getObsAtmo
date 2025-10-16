Rebuilding parameter grid
=========================


The atmospheric parameters - transmission grid may be regenerated
by using the python script ``rebuildGrid.py``.

Theses grids consist in a set of 5 files 1 set per pre-defined observation site.

* summary file containing data grid definitions under the format of a dictionnary(*pickle format*)
* the 2D (wavelength,airmass) grid for the Rayleigh scattering transmission,
* the 2D (wavelength,airmass) grid for the O2 absorption transmission,
* the 3D (wavelengthn,airmass, pwv) grid for the precipitable water vapor tranmission,
* the 3D (wavelengthn,airmass, ozone) grid for the Ozone transmission.

Requirements for rebuildGrid
----------------------------

Runnning ``rebuildGrid.py`` requires the installation of ``libRadtran``
and ``libradtranpy``.

Documentation to install those packages can be found at:

* `readthedocs site for libradtranpy  <https://libradtranpy.readthedocs.io/en/latest/>`_
* `github repository for libradtranpy <https://github.com/LSSTDESC/libradtranpy/>`_


Calling rebuildGrid
-------------------


In a shell :

| >> python rebuildGrids.py -h
| ************************************************************************
| Help to generation of atmospheric parameter grid for getObsAtmo emulator
| rebuildGrids.py  -s<observation site-string> -a <airmassmin,airmassmax,nbins> -v <pwvmin, pwvmax,nbins> -o <ozmin, ozmax,nbins>
| Observation sites are :
| LSST CTIO OHP PDM OMK OSL
| - atmospheric parameters should be specified of 3 numbers : valmin, valmax, x
|  where valmin is the minimum value, valmax is the maximum value
|        * if x = N in integer  -->  number of grid points
|        * if x = dx is float   -->  point spacing
|        -a  : airmass
|        -v  : precipitable water vapor
|        -o  : ozone
|        --------------------------------------------------------
| example : python  rebuildGrids.py  -s LSST -a 1,2.6,0.1 -v 0,15.25,0.25 -o 0.,625.0,25.
| or      : python  rebuildGrids.py  -s LSST -a 1,2.5,16   -v 0,15.25,61  -o  0.,625.0,25
|	 Arguments actually provided :
|	 	 Number of arguments: 2 arguments.
|	 	 Argument List: ['rebuildGrids.py', '-h']
