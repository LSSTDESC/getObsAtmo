# getObsAtmo

Python package to emulate atmospheric transparency simulation for different observation sites.

Transmission profiles depending on wavelength and airmass for Rayleigh scattering and atmospheric components absorption like Oxygen,
water vapor or Ozon were extracted from [libradtran](http://www.libradtran.org).
Those transmission profiles are located in `obsatmo_data`.
In addition to libradtran interpolated gridded profiles, an analytic expression scattering for single component aerosol scattering is provided.


# Download


```bash
git clone https://github.com/LSSTDESC/getObsAtmo.git
```


# Installation


``` bash
cd getObsAtmo
python setup.py install
```


# tests


```bash
python -m unittest
```

# Example


```python

from getObsAtmo.getObsAtmo import ObsAtmo

emul = ObsAtmo(obs_str = "LSST", pressure = 743.0)
wl = [400.,800.,900.]
am=1.2
pwv =4.0
oz=300.
transm = emul.GetAllTransparencies(wl,am,pwv,oz) # get transmission
print("wavelengths (nm) \t = ",wl)
print("transmissions    \t = ",transm)

emul.plot_transmission()  # plot the transmission
```



