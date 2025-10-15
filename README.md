# getObsAtmo

Python package to emulate atmospheric transparency simulation for different observation sites (version 0.2.0, October 2025).

Transmission profiles depending on wavelength and airmass for Rayleigh scattering and atmospheric components absorption like Oxygen,
water vapor or Ozon were extracted from [libradtran](http://www.libradtran.org).
Those transmission profiles are located in `getObsAtmo/obsatmo_data`.
In addition to libradtran interpolated gridded profiles, an analytic expression scattering for single component aerosol scattering is provided.

## Download

```bash
git clone https://github.com/LSSTDESC/getObsAtmo.git
```

## Installation

```bash
cd getObsAtmo
pip install .
```

or

```bash
pip install -e .
```

or if you want the documentation

```bash
pip install .[docs]
```

or

```bash
pip install -e .[docs]
```

## tests

```bash
python -m unittest
```

or

```bash
 python -m unittest tests/test_getObsAtmo.py
```

## Documentation

After the installation of sphinx aand `pip install .[docs]`, you can build the documentation with:

```bash
cd docs
make html
```

The documentation will be in `docs/build/html/index.html`

## Usage

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

## Example notebooks

Example notebooks are provided in the `examples/notebooks` directory.

## Online Documentation

The detailed documentation can be used in https://getobsatmo.readthedocs.io/en/latest/index.html
