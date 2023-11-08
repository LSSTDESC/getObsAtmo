#!/usr/bin/env python
# coding: utf-8

# # Generate grids


# Import some generally useful packages

import os
import numpy as np
from itertools import cycle, islice
import copy
import pickle

from scipy import interpolate
from libradtranpy import  libsimulateVisible
import sys,getopt


import warnings
warnings.filterwarnings('ignore')


# preselected observation sites 
Dict_Of_sitesAltitudes = {'LSST':2.663, # Rubin-LSST
                          'CTIO':2.207, # Cerro Tololo Inter-American Observatory
                          'OHP':0.65, # Observatoire de Haute Provence
                          'PDM':2.8905, # Observatoire du Pic du Midi
                          'OMK':4.205,  # Mauna Kea
                          'OSL':0,      # Sea Level
                           }

def usage():
    print("*******************************************************************")
    print(sys.argv[0],' -s<observation site-string> -a <airmassmin,airmassmax,nbins> -v <pwvmin, pwvmax,nbins> -o <ozmin, ozmax,nbins>')
    print("Observation sites are : ")
    print(' '.join(Dict_Of_sitesAltitudes.keys()))
    print('example : python ',  sys.argv[0], ' -s LSST -a 0,2.5,10 -v 0,10,20 -o 0,600,30' )
   
    print('\t Arguments actually provided : ')
    print('\t \t Number of arguments:', len(sys.argv), 'arguments.')
    print('\t \t Argument List:', str(sys.argv))

def decode_args(arg_str):
    return arg_str.split(',')

if __name__ == "__main__":


    try:
        opts, args = getopt.getopt(sys.argv[1:],"hs:a:v:o:",["s=","a=","v=","o="])
    except getopt.GetoptError:
        print(' Exception bad getopt with :: '+sys.argv[0]+ ' -s<observation-site-string> -a <airmassmin,airmassmax,nbins> -v <pwvmin,pwvmax,nbins> -o <ozmin,ozmax,nbins>' )
        sys.exit(2)

    print('opts = ',opts)
    print('args = ',args)
        
     
    alt_str = ""
    am_str  = ""
    pwv_str =""
    oz_str  = ""

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-s", "--site"):
            alt_str = arg.upper()
        elif opt in ("-a", "--airmass"):
            am_str = arg
        elif opt in ("-v", "--pwv"):
            pwv_str = arg
        elif opt in ("-o", "--oz"):
            oz_str = arg
        else:
            msg =f"Unknown option {opt} = {arg}"
        
    if alt_str in Dict_Of_sitesAltitudes.keys():
        OBS_tag = alt_str
    else:
        print(f"Observatory {alt_str} not in preselected observation site")
        print(f"This site {alt_str} must be added in libradtranpy preselected sites")
        sys.exit()


    print("am_str=",am_str,decode_args(am_str))


    print("pwv_str=",pwv_str,decode_args(pwv_str))
    print("oz_str=",oz_str,decode_args(oz_str))
   
    exit(0)
    
    FLAG_SCATTERING = True
    FLAG_ABSORPTION = True
   

    ########################
    ## Configuration
    ########################


    file0_out = f"{OBS_tag}_atmospherictransparencygrid_params.pickle"
    file1_out = f"{OBS_tag}_atmospherictransparencygrid_rayleigh.npy"
    file2_out = f"{OBS_tag}_atmospherictransparencygrid_O2abs.npy"
    file3_out = f"{OBS_tag}_atmospherictransparencygrid_PWVabs.npy"
    file4_out = f"{OBS_tag}_atmospherictransparencygrid_OZabs.npy"
    

    ####################
    # info dictionary

    info_params = {}


    info_params["OBS"] = OBS_tag

    ##################
    # ### wavelength
    ##################

    WLMIN=300.
    WLMAX=1100.
    WLBIN=1.
    NWLBIN=int((WLMAX-WLMIN)/WLBIN)
    WL=np.linspace(WLMIN,WLMAX,NWLBIN)

    info_params["WLMIN"] = WLMIN
    info_params["WLMAX"] = WLMAX
    info_params["WLBIN"] = WLBIN
    info_params["NWLBIN"] = NWLBIN
    info_params["WL"] = WL

    #################################
    # training and test dictionaries
    #################################

    info_params = copy.deepcopy(info_params)
    

    ####################
    # ### airmass
    ####################

    AIRMASSMIN=1.0
    AIRMASSMAX=2.6
    DAM = 0.1

    airmasses = np.arange(AIRMASSMIN,AIRMASSMAX,DAM)
    NAM=len(airmasses)


    NX=len(airmasses)
    NY=NWLBIN

    info_params["AIRMASSMIN"] = airmasses.min()
    info_params["AIRMASSMAX"] = airmasses.max()
    info_params["NAIRMASS"] = len(airmasses)
    info_params["DAIRMASS"] = np.median(np.diff(airmasses))
    info_params["AIRMASS"]  = airmasses


    #########################
    # ### PWV
    #########################

    PWVMIN = 0
    PWVMAX = 11
    DPWV = 0.25

    pwvs = np.arange(PWVMIN,PWVMAX,DPWV)
    

    NPWV = len(pwvs)

    info_params["PWVMIN"] = pwvs.min()
    info_params["PWVMAX"] = pwvs.max()
    info_params["NPWV"] = len(pwvs)
    info_params["DPWV"] = np.median(np.diff(pwvs))
    info_params["PWV"]  = pwvs


    # ### OZONE
    OZMIN = 0
    OZMAX = 600
    DOZ   = 25

    ozs = np.arange(OZMIN,OZMAX,DOZ)
    
    NOZ = len(ozs)

    info_params["OZMIN"] = ozs.min()
    info_params["OZMAX"] = ozs.max()
    info_params["NOZ"] = len(ozs)
    info_params["DOZ"] = np.median(np.diff(ozs))
    info_params["OZ"]  = ozs

    
    with open(file0_out, 'wb') as handle:
        pickle.dump(info_params, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    ##############################
    # ### Data container initialisation
    ####################################

    data_O2abs = np.zeros((NWLBIN,NAM))
    data_rayleigh = np.zeros((NWLBIN,NAM))
    data_H2Oabs = np.zeros((NWLBIN,NAM,NPWV))
    data_OZabs = np.zeros((NWLBIN,NAM,NOZ))
    

    ##########################################
    # Simulation of Rayleigh scattering
    ##############################################

    # note we call function ProcessSimulation with proc_str='sc' which mean scattering only

    # set pwv and ozone to zero (no need here)

    if FLAG_SCATTERING:
        pwv= 0
        oz = 0

        #training
        for idx,am in enumerate(airmasses):
            path,thefile = libsimulateVisible.ProcessSimulation(am,pwv,oz,0,prof_str='us',proc_str='sc',cloudext=0.0, altitude_str = OBS_tag,FLAG_VERBOSE=False)
            data = np.loadtxt(os.path.join(path,thefile))
            f = interpolate.interp1d(x=data[:,0], y=data[:,1],fill_value="extrapolate")
            atm=f(WL)
            data_rayleigh[:,idx]=atm
    
        np.save(file1_out,data_rayleigh, allow_pickle=False)



    if FLAG_ABSORPTION:
        ##########################################
        # Simulation of O2 absorption
        ##############################################

        # note we call function ProcessSimulation with proc_str='ab' which mean absorption only

        # set pwv and ozone to zero (no need here)
        pwv = 0.0
        oz = 0.0


        print("======================================")
        print("Simulation of O2 sample")
        print("======================================")

        for idx,am in enumerate(airmasses):
            path,thefile = libsimulateVisible.ProcessSimulation(am,pwv,oz,0,prof_str='us',proc_str='ab',cloudext=0.0,altitude_str = OBS_tag ,FLAG_VERBOSE=False)
            data = np.loadtxt(os.path.join(path,thefile))
            f = interpolate.interp1d(x=data[:,0], y=data[:,1],fill_value="extrapolate")
            atm=f(WL)
            data_O2abs[:,idx]=atm

        np.save(file2_out,data_O2abs, allow_pickle=False)
        print(f"...... O2 abs file {file2_out} written")

        

        ##########################################
        # Simulation of H2O absorption
        ##############################################

        # ## Precipitable water vapor
        print("======================================")
        print("Simulation of PWV sample")
        print("======================================")

        # note we call function ProcessSimulation with proc_str='ab' which mean absorption only
        # we cut ozone
        # however we cannot cut O2 absorption

        oz = 0.0
        for idx_pwv,pwv in enumerate(pwvs):
            data_slice=np.zeros((NWLBIN,NAM))
            for idx_am,am in enumerate(airmasses):     
                path,thefile = libsimulateVisible.ProcessSimulation(am,pwv,oz,0,prof_str='us',proc_str='ab',cloudext=0.0, altitude_str = OBS_tag ,FLAG_VERBOSE=False)
                data = np.loadtxt(os.path.join(path,thefile))
                f = interpolate.interp1d(x=data[:,0], y=data[:,1],fill_value="extrapolate")
                atm=f(WL)
                data_slice[:,idx_am]=atm
            # remove O2 absorption in H2O absorption profile
            data_slice /=data_O2abs
            data_H2Oabs[:,:,idx_pwv] = data_slice

        np.save(file3_out,data_H2Oabs,allow_pickle=False)
        print(f"...... H2O abs file {file3_out} written")
       


        ##########################################
        # Simulation of Ozone absorption
        ##############################################

        print("======================================")
        print("Simulation of Ozone sample")
        print("======================================")

        # note we call function ProcessSimulation with proc_str='ab' which mean absorption only
        # we remove pwv absorption
        # however we cannot cut O2 absorption
        pwv=0.0
        for idx_oz,oz in enumerate(ozs):
            data_slice=np.zeros((NWLBIN,NAM))
            for idx_am,am in enumerate(airmasses):     
                path,thefile = libsimulateVisible.ProcessSimulation(am,pwv,oz,0,prof_str='us',proc_str='ab',cloudext=0.0, altitude_str = OBS_tag ,FLAG_VERBOSE=False)
                data = np.loadtxt(os.path.join(path,thefile))
                f = interpolate.interp1d(x=data[:,0], y=data[:,1],fill_value="extrapolate")
                atm=f(WL)
                data_slice[:,idx_am]=atm
            # remove O2 profile from ozone profile
            data_slice/=data_O2abs
            data_OZabs[:,:,idx_oz] = data_slice


        np.save(file4_out,data_OZabs, allow_pickle=False)
        print(f"...... O3 abs file {file4_out} written")

        
