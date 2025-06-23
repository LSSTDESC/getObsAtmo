# README.md

# getObsAtmo

This package provides a function functions to simulate atmospheric transmissions for different observatories in visible, near-infrared by sampling transmissions from libradtran.

## rebuildGrids:

- requires libradtran installed as well as libradtranpy installed.

- generate atmospheric transmissions grids required by getobsatmo
- usage: `rebuildGrids.py -s <site> -a <aerosol> -v <vis> -o <nir>`
- `-s` specifies the site, e.g., `LSST`, `OHP`, `ZTF`, `VLT`, `PDM`, `OSL`, `CTIO`, `OMK`
- `-a` specifies the aerosol model, e.g., `1,2.7,0.1` (default)
- `-v` specifies the visible wavelength range, e.g., `0,20.50,0.25` (default)
- `-o` specifies the near-infrared wavelength range, e.g., `0.,650.0,25.` (default)

````bash

```code
python  rebuildGrids.py  -s LSST -a 1,2.7,0.1 -v 0,20.50,0.25 -o 0.,650.0,25.
python  rebuildGrids.py  -s OHP -a 1,2.7,0.1 -v 0,25.50,0.25 -o 0.,650.0,25.
python  rebuildGrids.py  -s ZTF -a 1,2.7,0.1 -v 0,30.50,0.25 -o 0.,650.0,25.
python  rebuildGrids.py  -s VLT -a 1,2.7,0.1 -v 0,25.50,0.25 -o 0.,650.0,25.
python  rebuildGrids.py  -s PDM -a 1,2.7,0.1 -v 0,20.50,0.25 -o 0.,650.0,25.
python  rebuildGrids.py  -s OSL -a 1,2.7,0.1 -v 0,25.50,0.25 -o 0.,650.0,25.
python  rebuildGrids.py  -s CTIO -a 1,2.7,0.1 -v 0,20.50,0.25 -o 0.,650.0,25.
python  rebuildGrids.py  -s OMK -a 1,2.7,0.1 -v 0,15.50,0.25 -o 0.,650.0,25.000

````
