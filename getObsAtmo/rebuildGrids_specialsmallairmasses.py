#!/usr/bin/env python
# coding: utf-8
# # Generate atmospheric parameter grids of transmissions 
# for getObsAtmo emulator
# Dedicated to sample finely small airmasses because of some artifact close to z=0
# last update 2024-11-19
####################################

# Import some generally useful packages

import os
import numpy as np
from itertools import cycle, islice
import copy
import pickle

from scipy import interpolate
from libradtranpy import  libsimulateVisible
import sys,getopt

import ast
import numbers 

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
    print("************************************************************************")
    print("Help to generation of atmospheric parameter grid for getObsAtmo emulator")
    print(sys.argv[0],' -s<observation site-string> -a <airmassmin,airmassmax,nbins> -v <pwvmin, pwvmax,nbins> -o <ozmin, ozmax,nbins>')
    print("Observation sites are : ")
    print(' '.join(Dict_Of_sitesAltitudes.keys()))
    print('- atmospheric parameters should be specified of 3 numbers : valmin, valmax, x' )
    print('  where valmin is the minimum value, valmax is the maximum value ')
    print('        * if x = N in integer  -->  number of grid points')
    print('        * if x = dx is float   -->  point spacing ')
    print('        -a  : airmass ')
    print('        -v  : precipitable water vapor ')
    print('        -o  : ozone ')
    print('        --------------------------------------------------------')
    print(' example : python ',  sys.argv[0], ' -s LSST -a 1,2.6,0.1 -v 0,15.25,0.25 -o 0.,625.0,25.' )
    print(' or      : python ',  sys.argv[0], ' -s LSST -a 1,2.5,16  -v 0,15.25,61  -o  0.,625.0,25 ' )
    print(" >>>> for specialsmallairmass, please don't provide the a parameter : a special grid is generated")
    print('\t Arguments actually provided : ')
    print('\t \t Number of arguments:', len(sys.argv), 'arguments.')
    print('\t \t Argument List:', str(sys.argv))
    print("************************************************************************")

def is_float(str):
    try:
        float(str)
        return True
    except ValueError:
        return False

def decode_args(arg_str):
    """String the input args in three numbers 
    :param arg_str: 
    :type arg_str: _type_
    :raises Exception: either the input arg cannot be splint in 3 element, or one of element is not a number
    :return: a numpy array of numbers
    :rtype: float
    """
    list_of_numbers = arg_str.split(',')

    nn = len(list_of_numbers)
    if nn != 3:
        msg = f"bad list of 3 numbers {arg_str}" 
        raise Exception(msg)
    
    test_all_floats = [ is_float(str) for str in list_of_numbers ]
    res_floats = all(test for test in test_all_floats) 
    test_all_digits = [ str.isdigit() for str in list_of_numbers ]    
    res_digits = all(test for test in test_all_digits) 

    if not res_floats:
        msg = f"bad list of 3 numbers {arg_str} are not all numbers" 
        raise Exception(msg)
    
    #print("test floats : ",  test_all_floats , "==>" , res_floats)
    #print("test digits : ",  test_all_digits , "==>" , res_digits)
    
    valmin= float(list_of_numbers[0])
    valmax= float(list_of_numbers[1])

    if valmin>= valmax:
        msg = f"bad values for {arg_str} : valmin = {valmin} >= valmax = {valmax}" 
        raise Exception(msg)
    
    if  test_all_digits[2]:
        N = int(list_of_numbers[2])
        array = np.linspace(valmin,valmax,N)
        return array
    else:
        step = float(list_of_numbers[2])
        array = np.arange(valmin,valmax,step)
        return array


    
#-----------------------------------------------------------------------
if __name__ == "__main__":


    try:
        opts, args = getopt.getopt(sys.argv[1:],"hs:a:v:o:",["s=","a=","v=","o="])
    except getopt.GetoptError:
        print(' Exception bad getopt with :: '+sys.argv[0]+ ' -s<observation-site-string> -a <airmassmin,airmassmax,nbins> -v <pwvmin,pwvmax,nbins> -o <ozmin,ozmax,nbins>' )
        sys.exit(2)

    #print('opts = ',opts)
    #print('args = ',args)
        
     
    alt_str = ""
    am_str  = ""
    pwv_str = ""
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

    if am_str != "":
        am_array = decode_args(am_str)
    else:
        am_array=np.array([])

    if pwv_str != "":
        pwv_array = decode_args(pwv_str)
    else:
        pwv_array=np.array([])
    
    if oz_str != "":
        oz_array = decode_args(oz_str)
    else:
        oz_array=np.array([])

    

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
    #  dictionary
    #################################

    info_params = copy.deepcopy(info_params)
    
    ####################
    # ### airmass
    ####################

    print("-------------------------------------------------------------------------------------------------")

    if am_str == "":
        AIRMASSMIN=1.0
        AIRMASSMAX=2.51
        DAM = 0.1
        airmasses = np.arange(AIRMASSMIN,AIRMASSMAX,DAM)
        airmass_small = np.arange(1,1.2,0.01)
        airmasses=np.union1d(airmass_small,airmasses)
        AIRMASSMAX = airmasses.max()
        print(">>> Select default airmass-grid : ")
    else:
        AIRMASSMIN = am_array.min()
        AIRMASSMAX = am_array.max()
        DAM = np.median(np.diff(am_array))
        airmasses = am_array
        print(">>> Select user airmass-grid : ")

    NAM=len(airmasses)
    print(f"NAM = {NAM} : ",airmasses)


    #NX=len(airmasses)
    #NY=NWLBIN

    info_params["AIRMASSMIN"] = airmasses.min()
    info_params["AIRMASSMAX"] = airmasses.max()
    info_params["NAIRMASS"] = len(airmasses)
    info_params["DAIRMASS"] = np.median(np.diff(airmasses))
    info_params["AIRMASS"]  = airmasses

    msg = f"airmass - grid : Npoints = {NAM}, valmin={AIRMASSMIN:.3f} , valmax={AIRMASSMAX:.3f} , valstep= {DAM:.3f}"
    print(msg)
    print(airmasses)

    #########################
    # ### PWV
    #########################

    print("-------------------------------------------------------------------------------------------------")

    if pwv_str == "":
        PWVMIN = 0.0
        PWVMAX = 25.50
        DPWV = 0.25
        pwvs = np.arange(PWVMIN,PWVMAX,DPWV)
        PWVMAX =pwvs.max()
        print(">>> Select default pwv-grid : ")
    else:
        PWVMIN = pwv_array.min()
        PWVMAX = pwv_array.max()
        DPWV = np.median(np.diff(pwv_array))
        pwvs = pwv_array
        print(">>> Select user pwv-grid : ")


    NPWV = len(pwvs)

    info_params["PWVMIN"] = pwvs.min()
    info_params["PWVMAX"] = pwvs.max()
    info_params["NPWV"] = len(pwvs)
    info_params["DPWV"] = np.median(np.diff(pwvs))
    info_params["PWV"]  = pwvs

    msg = f"pwv - grid : Npoints = {NPWV}, valmin={PWVMIN:.3f} , valmax={PWVMAX:.3f} , valstep= {DPWV:.3f}"
    print(msg)
    print(pwvs)


    # ### OZONE
    print("-------------------------------------------------------------------------------------------------")

    if oz_str == "":
        OZMIN = 0.0
        OZMAX = 610.0
        DOZ   = 25.0
        ozs = np.arange(OZMIN,OZMAX,DOZ)
        OZMAX = ozs.max()
        print(">>> Select default oz-grid : ")
    else:
        OZMIN = oz_array.min()
        OZMAX = oz_array.max()
        DOZ = np.median(np.diff(oz_array))
        ozs = oz_array
        print(">>> Select user oz-grid : ")

    
    NOZ = len(ozs)

    info_params["OZMIN"] = ozs.min()
    info_params["OZMAX"] = ozs.max()
    info_params["NOZ"] = len(ozs)
    info_params["DOZ"] = np.median(np.diff(ozs))
    info_params["OZ"]  = ozs

    msg = f"oz - grid : Npoints = {NOZ}, valmin={OZMIN:.3f} , valmax={OZMAX:.3f} , valstep= {DOZ:.3f}"
    print(msg)
    print(ozs)

    print("-------------------------------------------------------------------------------------------------")


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
    pwv= 0
    oz = 0

    
   
    for idx,am in enumerate(airmasses):
        wl, atm = libsimulateVisible.ProcessSimulation(am,pwv,oz,press_num=0, aer_num=0,angstrom_exponent_num=1.4,prof_str='us',proc_str='sc',cloudext=0.0, altitude = OBS_tag,aer_lambda0=500.,FLAG_VERBOSE=False)
        f = interpolate.interp1d(x=wl, y=atm,fill_value="extrapolate")
        atm=f(WL)
        data_rayleigh[:,idx]=atm
    
    np.save(file1_out,data_rayleigh, allow_pickle=False)

    
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
        wl, atm = libsimulateVisible.ProcessSimulation(am,pwv,oz,press_num=0, aer_num=0,angstrom_exponent_num=1.4,prof_str='us',proc_str='ab',cloudext=0.0,altitude = OBS_tag ,FLAG_VERBOSE=False)
        f = interpolate.interp1d(x=wl, y=atm,fill_value="extrapolate")
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
            wl,atm = libsimulateVisible.ProcessSimulation(am,pwv,oz,press_num=0, aer_num=0,angstrom_exponent_num=1.4,prof_str='us',proc_str='ab',cloudext=0.0, altitude = OBS_tag ,FLAG_VERBOSE=False)
            f = interpolate.interp1d(x=wl, y=atm,fill_value="extrapolate")
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
    print(" Simulation of Ozone sample" )
    print("======================================")

    # note we call function ProcessSimulation with proc_str='ab' which mean absorption only
    # we remove pwv absorption
    # however we cannot cut O2 absorption
    pwv=0.0
    for idx_oz,oz in enumerate(ozs):
        data_slice=np.zeros((NWLBIN,NAM))
        for idx_am,am in enumerate(airmasses):     
            wl,atm = libsimulateVisible.ProcessSimulation(am,pwv,oz,press_num=0, aer_num=0,angstrom_exponent_num=1.4,prof_str='us',proc_str='ab',cloudext=0.0, altitude = OBS_tag ,FLAG_VERBOSE=False)
            f = interpolate.interp1d(x=wl, y=atm,fill_value="extrapolate")
            atm=f(WL)
            data_slice[:,idx_am]=atm
        # remove O2 profile from ozone profile
        data_slice/=data_O2abs
        data_OZabs[:,:,idx_oz] = data_slice

    np.save(file4_out,data_OZabs, allow_pickle=False)
    print(f"...... O3 abs file {file4_out} written")

        
