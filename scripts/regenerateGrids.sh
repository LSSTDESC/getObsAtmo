  #! bin/sh
# This script regenerates the grids for the getObsAtmo package.
# It should be run from the root of the getObsAtmo package directory.
# Usage: ./getObsAtmo/scripts/regenerateGrids.sh


#python ../getObsAtmo/rebuildGrids.py -s LSST -a 1,2.7,0.1 -v 0,20.50,0.25 -o 0.,650.0,25.
#python ../getObsAtmo/rebuildGrids.py -s OHP -a 1,2.7,0.1 -v 0,25.50,0.25 -o 0.,650.0,25.
python ../getObsAtmo/rebuildGrids.py -s ZTF -a 1,2.7,0.1 -v 0,30.50,0.25 -o 0.,650.0,25.
#python ../getObsAtmo/rebuildGrids.py -s VLT -a 1,2.7,0.1 -v 0,25.50,0.25 -o 0.,650.0,25.
#python ../getObsAtmo/rebuildGrids.py -s PDM -a 1,2.7,0.1 -v 0,20.50,0.25 -o 0.,650.0,25.
#python ../getObsAtmo/rebuildGrids.py -s OSL -a 1,2.7,0.1 -v 0,25.50,0.25 -o 0.,650.0,25.
#python ../getObsAtmo/rebuildGrids.py -s CTIO -a 1,2.7,0.1 -v 0,20.50,0.25 -o 0.,650.0,25.
#python ../getObsAtmo/rebuildGrids.py -s OMK -a 1,2.7,0.1 -v 0,15.50,0.25 -o 0.,650.0,25.000