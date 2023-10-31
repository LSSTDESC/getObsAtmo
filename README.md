# getObsAtmo
Python package to emulate atmospheric transparency simulation for different observation sites.

Transmission profiles depending on wavelength and airmass for Rayleigh scattering and atmospheric components absorption like Oxygen,
water vapor or Ozon were extracted from [libradtran](http://www.libradtran.org).
Those transmission profiles are located in `obsatmo_data`.
In addition to libradtran profiles, a analytic scattering profile aerosol for a single component is provided.



Example:
```
from getObsAtmo.getObsAtmo import *

test = is_calspec("eta1 dor")
c = Calspec("eta1 dor")
c.get_spectrum_fits_filename()  # download the fits file
c.get_spectrum_table()  # download and return an Astropy table
c.get_spectrum_numpy()  # download and return a dictionnary of numpy arrays with units
c.plot_spectrum()  # download and plot the spectrum
```
