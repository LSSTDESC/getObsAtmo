#!/usr/bin/env python
# coding: utf-8
# # Generate atmospheric parameter grids of transmissions
# for getObsAtmo emulator
####################################

# Import some generally useful packages

import os
import numpy as np
from itertools import cycle, islice
import copy
import pickle
import argparse

from scipy import interpolate
from libradtranpy import libsimulateVisible
import sys, getopt
from joblib import Parallel, delayed

import ast
import numbers

import warnings
import time

warnings.filterwarnings("ignore")


# preselected observation sites
Dict_Of_sitesAltitudes = {
    "LSST": 2.663,  # Rubin-LSST
    "CTIO": 2.207,  # Cerro Tololo Inter-American Observatory
    "OHP": 0.65,  # Observatoire de Haute Provence
    "PDM": 2.8905,  # Observatoire du Pic du Midi
    "OMK": 4.205,  # Mauna Kea
    "OSL": 0,  # Sea Level
}


# pressure calculated by libradtran
Dict_Of_sitesPressures = {
    "LSST": 731.50433,
    "CTIO": 774.6052,
    "OHP": 937.22595,
    "PDM": 710.90637,
    "OMK": 600.17224,
    "OSL": 1013.000,
}


def usage():
    print("************************************************************************")
    print("Help to generation of atmospheric parameter grid for getObsAtmo emulator")
    print(
        sys.argv[0],
        " -s<observation site-string> -a <airmassmin,airmassmax,nbins> -v <pwvmin, pwvmax,nbins> -o <ozmin, ozmax,nbins>",
    )
    print("Observation sites are : ")
    print(" ".join(Dict_Of_sitesAltitudes.keys()))
    print(
        "- atmospheric parameters should be specified of 3 numbers : valmin, valmax, x"
    )
    print("  where valmin is the minimum value, valmax is the maximum value ")
    print("        * if x = N in integer  -->  number of grid points")
    print("        * if x = dx is float   -->  point spacing ")
    print("        -a  : airmass ")
    print("        -v  : precipitable water vapor ")
    print("        -o  : ozone ")
    print("        --------------------------------------------------------")
    print(
        " example : python ",
        sys.argv[0],
        " -s LSST -a 1,2.6,0.1 -v 0,10.25,0.25 -o 0.,625.0,25.",
    )
    print(
        " or      : python ",
        sys.argv[0],
        " -s LSST -a 1,2.5,16   -v 0,10,41  -o  0.,600.0,25 ",
    )
    print("\t Arguments actually provided : ")
    print("\t \t Number of arguments:", len(sys.argv), "arguments.")
    print("\t \t Argument List:", str(sys.argv))
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
    list_of_numbers = arg_str.split(",")

    nn = len(list_of_numbers)
    if nn != 3:
        msg = f"bad list of 3 numbers {arg_str}"
        raise Exception(msg)

    test_all_floats = [is_float(str) for str in list_of_numbers]
    res_floats = all(test for test in test_all_floats)
    test_all_digits = [str.isdigit() for str in list_of_numbers]
    res_digits = all(test for test in test_all_digits)

    if not res_floats:
        msg = f"bad list of 3 numbers {arg_str} are not all numbers"
        raise Exception(msg)

    # print("test floats : ",  test_all_floats , "==>" , res_floats)
    # print("test digits : ",  test_all_digits , "==>" , res_digits)

    valmin = float(list_of_numbers[0])
    valmax = float(list_of_numbers[1])

    if valmin >= valmax:
        msg = f"bad values for {arg_str} : valmin = {valmin} >= valmax = {valmax}"
        raise Exception(msg)

    if test_all_digits[2]:
        N = int(list_of_numbers[2])
        array = np.linspace(valmin, valmax, N)
        return array
    else:
        step = float(list_of_numbers[2])
        array = np.arange(valmin, valmax, step)
        return array


# -----------------------------------------------------------------------
if __name__ == "__main__":

    # Create the parser
    parser = argparse.ArgumentParser(description="Process some parameters.")

    # Define the arguments
    parser.add_argument("-s", "--site", type=str, help="Observation site string")
    parser.add_argument("-a", "--airmass", type=str, help="Airmass min, max, nbins")
    parser.add_argument("-v", "--pwv", type=str, help="PWV min, max, nbins")
    parser.add_argument("-o", "--oz", type=str, help="Ozone min, max, nbins")
    parser.add_argument(
        "-p", "--profile", default="afglus", type=str, help="Profile name"
    )
    parser.add_argument(
        "-r",
        "--resolution",
        default="coarse",
        required=False,
        type=str,
        help="Resolution of the simulation",
    )
    parser.add_argument(
        "-wlmin", "--wlmin", default=300, required=False, type=int, help="Min lambda"
    )
    parser.add_argument(
        "-wlmax", "--wlmax", default=1200, required=False, type=int, help="Max lambda"
    )
    parser.add_argument(
        "-wlbin",
        "--wlbin",
        default=0.5,
        required=False,
        type=float,
        help="Lambda bin for the interpolated spectra",
    )
    parser.add_argument(
        "-rs",
        "--rte",
        default="twostr",
        required=False,
        type=str,
        help="Radiative transfer equation solver",
    )

    # Parse the arguments
    args = parser.parse_args()

    # Assign the arguments to variables
    alt_str = args.site.upper() if args.site else ""
    am_str = args.airmass if args.airmass else ""
    pwv_str = args.pwv if args.pwv else ""
    oz_str = args.oz if args.oz else ""
    atm_model = args.profile
    molresol = args.resolution
    wlmin = args.wlmin
    wlmax = args.wlmax
    wlbin = args.wlbin
    rte = args.rte

    t_start = time.time()

    # Check altitude str
    if alt_str in Dict_Of_sitesAltitudes.keys():
        OBS_tag = alt_str
        press_num = Dict_Of_sitesPressures[alt_str]
    else:
        print(f"Observatory {alt_str} not in preselected observation site")
        print(f"This site {alt_str} must be added in libradtranpy preselected sites")
        sys.exit()

    # Decode arguments for (min, max, step) of given parameters
    if am_str != "":
        am_array = decode_args(am_str)
    else:
        am_array = np.array([])

    if pwv_str != "":
        pwv_array = decode_args(pwv_str)
    else:
        pwv_array = np.array([])

    if oz_str != "":
        oz_array = decode_args(oz_str)
    else:
        oz_array = np.array([])

    ########################
    ## Configuration
    ########################

    # Filenames
    file0_out = (
        f"{OBS_tag}_{atm_model}_atmospherictransparencygrid_params.pickle"
    )
    file1_out = (
        f"{OBS_tag}_{atm_model}_atmospherictransparencygrid_rayleigh.npy"
    )
    file2_out = (
        f"{OBS_tag}_{atm_model}_atmospherictransparencygrid_O2abs.npy"
    )
    file3_out = (
        f"{OBS_tag}_{atm_model}_atmospherictransparencygrid_PWVabs.npy"
    )
    file4_out = (
        f"{OBS_tag}_{atm_model}_atmospherictransparencygrid_OZabs.npy"
    )

    ####################
    # info dictionary

    info_params = {}
    info_params["OBS"] = OBS_tag

    ##################
    # ### wavelength
    ##################

    WLMIN = wlmin
    WLMAX = wlmax
    WLBIN = wlbin
    NWLBIN = int((WLMAX - WLMIN) / WLBIN)
    WL = np.linspace(WLMIN, WLMAX, NWLBIN)

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

    print(
        "-------------------------------------------------------------------------------------------------"
    )

    if am_str == "":
        AIRMASSMIN = 1.0
        AIRMASSMAX = 2.51
        DAM = 0.1
        airmasses = np.arange(AIRMASSMIN, AIRMASSMAX, DAM)
        AIRMASSMAX = airmasses.max()
        print(">>> Select default airmass-grid : ")
    else:
        AIRMASSMIN = am_array.min()
        AIRMASSMAX = am_array.max()
        DAM = np.median(np.diff(am_array))
        airmasses = am_array
        print(">>> Select user airmass-grid : ")

    NAM = len(airmasses)

    # NX=len(airmasses)
    # NY=NWLBIN

    info_params["AIRMASSMIN"] = airmasses.min()
    info_params["AIRMASSMAX"] = airmasses.max()
    info_params["NAIRMASS"] = len(airmasses)
    info_params["DAIRMASS"] = np.median(np.diff(airmasses))
    info_params["AIRMASS"] = airmasses

    msg = f"airmass - grid : Npoints = {NAM}, valmin={AIRMASSMIN:.3f} , valmax={AIRMASSMAX:.3f} , valstep= {DAM:.3f}"
    print(msg)
    print(airmasses)

    #########################
    # ### PWV
    #########################

    print(
        "-------------------------------------------------------------------------------------------------"
    )

    if pwv_str == "":
        PWVMIN = 0.0
        PWVMAX = 10.24
        DPWV = 0.25
        pwvs = np.arange(PWVMIN, PWVMAX, DPWV)
        PWVMAX = pwvs.max()
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
    info_params["PWV"] = pwvs

    msg = f"pwv - grid : Npoints = {NPWV}, valmin={PWVMIN:.3f} , valmax={PWVMAX:.3f} , valstep= {DPWV:.3f}"
    print(msg)
    print(pwvs)

    # ### OZONE
    print(
        "-------------------------------------------------------------------------------------------------"
    )

    if oz_str == "":
        OZMIN = 0.0
        OZMAX = 610.0
        DOZ = 25.0
        ozs = np.arange(OZMIN, OZMAX, DOZ)
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
    info_params["OZ"] = ozs

    msg = f"oz - grid : Npoints = {NOZ}, valmin={OZMIN:.3f} , valmax={OZMAX:.3f} , valstep= {DOZ:.3f}"
    print(msg)
    print(ozs)

    print(
        "-------------------------------------------------------------------------------------------------"
    )

    with open(file0_out, "wb") as handle:
        pickle.dump(info_params, handle, protocol=pickle.HIGHEST_PROTOCOL)

    ##############################
    # ### Data container initialisation
    ####################################

    data_O2abs = np.zeros((NWLBIN, NAM))
    data_rayleigh = np.zeros((NWLBIN, NAM))
    data_H2Oabs = np.zeros((NWLBIN, NAM, NPWV))
    data_OZabs = np.zeros((NWLBIN, NAM, NOZ))

    def process_simulation(airmass, pwv, ozone, pressure, atm_model, proc):
        wl, atm = libsimulateVisible.ProcessSimulation(
            airmass=airmass,
            pwv=pwv,
            ozone=ozone,
            pressure=pressure,
            atm_model=atm_model,
            proc=proc,
            rte=rte,
            altitude=OBS_tag,
            lambda_min=WLMIN,
            lambda_max=WLMAX,
            wl_res=molresol,
            pseudospherical=True,
        )
        # data = np.loadtxt(os.path.join(path, thefile))
        f = interpolate.interp1d(x=wl, y=atm, fill_value="extrapolate")
        atm = f(WL)
        return atm

    ##########################################
    # Simulation of Rayleigh scattering
    ##############################################

    pwv = 0
    oz = 0

    results_rayleigh = Parallel(n_jobs=-1, backend="multiprocessing")(
        delayed(process_simulation)(
            airmass=am,
            pwv=pwv,
            ozone=oz,
            pressure=press_num,
            atm_model=atm_model,
            proc="no_absorption",
        )
        for am in airmasses
    )
    for idx, atm in enumerate(results_rayleigh):
        print(f"ITERATION {idx+1} / {len(airmasses)} FOR RAYLEIGH")
        data_rayleigh[:, idx] = atm

    np.save(file1_out, data_rayleigh, allow_pickle=False)

    ##########################################
    # Simulation of O2 absorption
    ##############################################

    pwv = 0.0
    oz = 0.0

    print("======================================")
    print("Simulation of O2 sample")
    print("======================================")

    results_O2abs = Parallel(n_jobs=-1, backend="multiprocessing")(
        delayed(process_simulation)(
            airmass=am,
            pwv=pwv,
            ozone=oz,
            pressure=press_num,
            atm_model=atm_model,
            proc="no_scattering",
        )
        for am in airmasses
    )
    for idx, atm in enumerate(results_O2abs):
        print(f"ITERATION {idx+1} / {len(airmasses)} FOR O2")
        data_O2abs[:, idx] = atm

    np.save(file2_out, data_O2abs, allow_pickle=False)
    print(f"...... O2 abs file {file2_out} written")

    ##########################################
    # Simulation of H2O absorption
    ##############################################

    print("======================================")
    print("Simulation of PWV sample")
    print("======================================")

    oz = 0.0
    for idx_pwv, pwv in enumerate(pwvs):
        print(f"ITERATION {idx+1} / {len(pwvs)} FOR PWV")
        results_H2Oabs = Parallel(n_jobs=-1, backend="multiprocessing")(
            delayed(process_simulation)(
                airmass=am,
                pwv=pwv,
                ozone=oz,
                pressure=press_num,
                atm_model=atm_model,
                proc="no_scattering",
            )
            for am in airmasses
        )
        data_slice = np.array(results_H2Oabs).T
        data_slice /= data_O2abs  # remove O2 absorption in H2O absorption profile
        data_H2Oabs[:, :, idx_pwv] = data_slice

    np.save(file3_out, data_H2Oabs, allow_pickle=False)
    print(f"...... H2O abs file {file3_out} written")

    ##########################################
    # Simulation of Ozone absorption
    ##############################################

    print("======================================")
    print(" Simulation of Ozone sample")
    print("======================================")

    pwv = 0.0
    for idx_oz, oz in enumerate(ozs):
        print(f"ITERATION {idx+1} / {len(pwvs)} FOR OZONE")
        results_OZabs = Parallel(n_jobs=-1, backend="multiprocessing")(
            delayed(process_simulation)(
                airmass=am,
                pwv=pwv,
                ozone=oz,
                pressure=press_num,
                atm_model=atm_model,
                proc="no_scattering",
            )
            for am in airmasses
        )
        data_slice = np.array(results_OZabs).T
        data_slice /= data_O2abs  # remove O2 profile from ozone profile
        data_OZabs[:, :, idx_oz] = data_slice

    np.save(file4_out, data_OZabs, allow_pickle=False)
    print(f"...... O3 abs file {file4_out} written")
    t_end = time.time()
    print(f"TOTAL TIME TO GENERATE GRIDS = {t_end-t_start:.1f} seconds")
