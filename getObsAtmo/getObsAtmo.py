# last update 2025-10-15 : use json instead of pickle
import numpy as np
from scipy.interpolate import RegularGridInterpolator

# import pickle
import json
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
from typing import Any, Dict, Optional


__all__ = [
    "Dict_Of_sitesAltitudes",
    "Dict_Of_sitesPressures",
    "Dict_Of_sitesAliases",
    "file_data_dict",
    "_getPackageDir",
    "get_obssite_keys",
    "sanitizeString",
    "validateObsName",
    "is_obssite",
    "ObsAtmo",
    "ObsAtmoPressure",
    "ObsAtmoGrid",
]


# preselected sites
Dict_Of_sitesAltitudes = {
    "LSST": 2.663,  # Rubin-LSST
    "CTIO": 2.207,  # Cerro Tololo Inter-American Observatory
    "OHP": 0.65,  # Observatoire de Haute Provence
    "ZTF": 1.712,  # Palomar Observatory
    "VLT": 2.635,  # Cerro Paranal (ESO)
    "PDM": 2.8905,  # Observatoire du Pic du Midi
    "OMK": 4.205,  # Mauna Kea
    "OSL": 0,  # Sea Level
}
# pressure calculated by libradtran
Dict_Of_sitesPressures = {
    "LSST": 731.50433,
    "CTIO": 774.6052,
    "OHP": 937.22595,
    "ZTF": 823.60004,  # Palomar Observatory
    "VLT": 734.08038,  # Cerro Paranal (ESO)
    "PDM": 710.90637,
    "OMK": 600.17224,
    "OSL": 1013.000,
}

Dict_Of_sitesTags = {
    "LSST": "LS",
    "CTIO": "CT",
    "OHP": "OH",
    "ZTF": "ZT",
    "VLT": "VL",
    "PDM": "PM",
    "OMK": "MK",
    "OSL": "SL",
}

Dict_Of_sitesAliases = {
    "LSST": ["Rubin", "Rubin Observatory", "Auxtel"],
    "CTIO": ["Cerro Tololo"],
    "OHP": ["Observatoire de Haute Provence"],
    "PDM": ["Pic du Midi", "Observatoire du Pic du Midi"],
    "OMK": ["Mauna Kea", "Mauna Kea Observatory"],
    "OSL": ["Sea Level", "Sea Level Observatory"],
    "ZTF": ["Palomar", "Palomar Observatory"],
    "VLT": [
        "Cerro Paranal",
        "Very Large Telescope",
        "Paranal Observatory",
        "Paranal",
    ],
}


file_data_dict = {
    # "info": "atmospherictransparencygrid_params.pickle",
    "info": "atmospherictransparencygrid_params.json",
    "data_rayleigh": "atmospherictransparencygrid_rayleigh.npy",
    "data_o2abs": "atmospherictransparencygrid_O2abs.npy",
    "data_pwvabs": "atmospherictransparencygrid_PWVabs.npy",
    "data_ozabs": "atmospherictransparencygrid_OZabs.npy",
}


def _getPackageDir():
    """This method must live in the top level of this package, so if this
    moves to a utils file then the returned path will need to account for that.
    """
    dirname = os.path.dirname(__file__)
    return dirname


def convert_dict_to_json(
    data_dict: Dict[str, Any], data_json: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Convert a Python dictionary containing NumPy arrays into a
    JSON-serializable dictionary.

    This function iterates over all key-value pairs in the input dictionary.
    NumPy arrays are converted to lists, while scalar types
    (int, float, str, bool, None)
    are preserved as-is. The resulting dictionary can be safely serialized
    using `json.dumps`.

    Args:
        data_dict (Dict[str, Any]): Input dictionary, possibly containing
        NumPy arrays or scalar values.
        data_json (Optional[Dict[str, Any]]): Optional output dictionary.
        If provided, the converted data will be written into it.
        Otherwise, a new dictionary will be created.

    Returns:
        Dict[str, Any]: A dictionary where all NumPy arrays have been converted
        into lists, ready for JSON serialization.

    Raises:
        TypeError: If a value in `data_dict` is not a supported type
        (NumPy array, scalar, or None).

    Examples:
        >>> import numpy as np
        >>> data = {'array': np.array([1, 2, 3]), 'value': 42}
        >>> json_data = convert_dict_to_json(data)
        >>> print(json_data)
        {'array': [1, 2, 3], 'value': 42}
    """
    if not isinstance(data_dict, dict):
        raise TypeError(f"Expected a dictionary, got {type(data_dict).__name__}.")

    if data_json is None:
        data_json = {}

    for key, value in data_dict.items():
        if isinstance(value, np.ndarray):
            data_json[key] = value.tolist()
        elif isinstance(value, (int, float, str, bool)) or value is None:
            data_json[key] = value
        else:
            raise TypeError(
                f"Unsupported type for key '{key}': {type(value).__name__}. "
                "Only NumPy arrays and scalar types are supported."
            )

    return data_json


def convert_json_to_dict(data_json: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert a JSON-like dictionary into a Python dictionary where
    list values are converted to NumPy arrays.

    This function iterates through all key-value pairs in the input dictionary.
    If a value is a list, it will be converted to a NumPy array.
    All other values are preserved as-is.

    Args:
        data_json (Dict[str, Any]): Input JSON-compatible dictionary
        (e.g., parsed from a JSON file).

    Returns:
        Dict[str, Any]: A dictionary where list values have been converted
        to NumPy arrays.

    Raises:
        TypeError: If the input is not a dictionary.
    """
    if not isinstance(data_json, dict):
        raise TypeError(f"Expected a dictionary, got {type(data_json).__name__}.")

    data_dict: Dict[str, Any] = {}
    for key, value in data_json.items():
        if isinstance(value, list):
            try:
                data_dict[key] = np.array(value)
            except Exception as e:
                raise ValueError(
                    f"Failed to convert list at key '{key}' " f"to a NumPy array: {e}"
                )
        else:
            data_dict[key] = value

    return data_dict


def get_obssite_dataframe():
    """
    Provide the list of observatories which have transmission grids

    Example
    -------

    .. doctest::

        >>> print(getObsSiteDataFrame())
             altitude   pressure
        LSST    2.663  731.50433
        CTIO    2.207   774.6052
        OHP      0.65  937.22595
        PDM    2.8905  710.90637
        OMK     4.205  600.17224
        OSL       0.0     1013.0
    """
    df = pd.DataFrame(
        columns=["altitude", "pressure"], index=list(Dict_Of_sitesAltitudes.keys())
    )
    for key in Dict_Of_sitesAltitudes.keys():
        df.loc[key] = pd.Series(
            {
                "altitude": Dict_Of_sitesAltitudes[key],
                "pressure": Dict_Of_sitesPressures[key],
            }
        )
    return df


def sanitizeString(label) -> str:
    """This method sanitizes the site label."""
    return label.upper().replace(" ", "")


def validateObsName(obssitename) -> str | None:
    """Validate if the obsite name is a valid observation site

    :param obssitename: observatory site name including a possible alias
    :type obssitename: str
    :return: valid obssite key label
    :rtype: str or None
    """

    sitename = sanitizeString(obssitename)
    for key_site, listnames in Dict_Of_sitesAliases.items():
        sanit_listname = list(map(lambda x: sanitizeString(x), listnames))
        if sitename in sanit_listname or sitename == key_site:
            return key_site
    return None


def get_obssite_keys(obs_label):
    """Return the DataFrame keys if an observation site name corresponds
    to a an entry in the tables.

    :param obs_label:  The observation site name. Five sites have ben
    preselected.
    :type obs_label: string among 'LSST','CTIO','OHP','PDM','OMK','OSL'

    :return: the dataframe corresponding to the location of the 'obs_label'
    :rtype: pandas series

    Examples
    --------
    .. doctest::

        >>> get_obssite_keys("lsst")  #doctest: +ELLIPSIS
        0     True
        1    False...

    """
    label = sanitizeString(obs_label)
    df = get_obssite_dataframe()
    name_index = [name.upper() for name in df.index]
    if len(name_index) > 0:
        keys = pd.Series([False] * len(df))
        for idx, name in enumerate(name_index):
            keys[idx] = name == label
        return keys
    else:
        raise KeyError("No observation site label {obs_label} in config dictionaries")


def is_obssite(obs_label):
    """
    Test if an observation site label is among the list atmospheric grids.

    :param obs_label:  The observation site name. Five sites have been
    preselected.
    :type obs_label: string among 'LSST','CTIO','OHP','PDM','OMK','OSL'

    :return:  check if the observation site is in the predefined set of
    observation sites
    :rtype: bool
        True if the observation site is star has a table.

    Examples
    --------
    >>> is_obssite("LSST")
    True
    >>> is_obssite("OHP")
    True
    """

    obs_tag = validateObsName(obs_label)

    if obs_tag is None:
        return False
    else:
        return np.any(get_obssite_keys(obs_tag))


class ObsAtmoGrid:
    """ "
    Emulate Atmospheric Transparency above LSST from a data grids
    extracted from libradtran and analytical functions for aerosols.
    There are 3 grids:
    - 2D grid Rayleigh transmission vs (wavelength,airmass)
    - 2D grid O2 absorption vs  (wavelength,airmass)
    - 3D grid for PWV absorption vs (wavelength,airmass,PWV)
    - 3D grid for Ozone absorption vs (wavelength,airmass,Ozone)
    - Aerosol transmission for any number of components

    :param obs_label:  The observation site name. Five sites have been
    preselected.
    :type obs_label: string among 'LSST','CTIO','OHP','PDM','OMK','OSL'
    """

    def __init__(self, obs_str="LSST"):
        """
        Initialize the class for data point files from which the 2D and 3D grids
        are created.
        Interpolation are calculated from the scipy RegularGridInterpolator()
        function
        Both types of data : trainging data for normal interpolaton use and the
        test data used to check accuracy of the interpolation of data.

        :param obs_label:  The observation site name. Five sites have been
        preselected.
        :type obs_label: string among 'LSST','CTIO','OHP','PDM','OMK','OSL'

        :raise: the `obs_label` is not in the list of predefined
        observatory sites

        :return: the emulator obsject
        :rtype: object of class `ObsAtmoGrid`

        Examples
        --------
        >>> emulator = ObsAtmoGrid('LSST')

        """
        self.OBS_tag = ""

        obs_tag = validateObsName(obs_str)

        if obs_tag is None:
            raise ValueError(
                f"Observatory {obs_str} not in preselected "
                f"observation sites.\n "
                f"This site {obs_str} must be added in libradtranpy "
                f"preselected sites "
                f"and generate corresponding scattering and absorption "
                f"profiles."
            )
        else:
            # print(f"{obs_str} site name validated as {obs_tag} observatory")
            self.OBS_tag = obs_tag

        self.Name = (
            f"Atmospheric emulator ObsAtmoGrid for observation site " f" {self.OBS_tag}"
        )

        # construct the path of input data files
        # self.path = os.path.join(_getPackageDir(), '../obsatmo_data')
        self.path = os.path.join(_getPackageDir(), "obsatmo_data")
        self.info_params = {}

        # load all data files (training and test)
        filename = os.path.join(self.path, self.OBS_tag + "_" + file_data_dict["info"])

        # with open(filename, 'rb') as f:
        #    self.info_params = pickle.load(f)

        if not os.path.exists(filename):
            raise FileNotFoundError(f"Missing JSON file: {filename}")

        with open(filename, "r") as f:
            loaded_data_json = json.load(f)
            self.info_params = convert_json_to_dict(loaded_data_json)

        data_rayleigh = np.load(
            os.path.join(
                self.path, self.OBS_tag + "_" + file_data_dict["data_rayleigh"]
            )
        )
        data_O2abs = np.load(
            os.path.join(self.path, self.OBS_tag + "_" + file_data_dict["data_o2abs"])
        )
        data_PWVabs = np.load(
            os.path.join(self.path, self.OBS_tag + "_" + file_data_dict["data_pwvabs"])
        )
        data_OZabs = np.load(
            os.path.join(self.path, self.OBS_tag + "_" + file_data_dict["data_ozabs"])
        )

        # setup training dataset (those used for interpolation)
        self.WLMIN = self.info_params["WLMIN"]
        self.WLMAX = self.info_params["WLMAX"]
        self.WLBIN = self.info_params["WLBIN"]
        self.NWLBIN = self.info_params["NWLBIN"]
        self.WL = self.info_params["WL"]
        self.OBS = self.info_params["OBS"]

        self.AIRMASSMIN = self.info_params["AIRMASSMIN"]
        self.AIRMASSMAX = self.info_params["AIRMASSMAX"]
        self.NAIRMASS = self.info_params["NAIRMASS"]
        self.DAIRMASS = self.info_params["DAIRMASS"]
        self.AIRMASS = self.info_params["AIRMASS"]

        self.PWVMIN = self.info_params["PWVMIN"]
        self.PWVMAX = self.info_params["PWVMAX"]
        self.NPWV = self.info_params["NPWV"]
        self.DPWV = self.info_params["DPWV"]
        self.PWV = self.info_params["PWV"]

        self.OZMIN = self.info_params["OZMIN"]
        self.OZMAX = self.info_params["OZMAX"]
        self.NOZ = self.info_params["NOZ"]
        self.DOZ = self.info_params["DOZ"]
        self.OZ = self.info_params["OZ"]

        # constant parameters defined for aerosol formula
        self.lambda0 = 500.0
        self.tau0 = 1.0

        # interpolation is done over training dataset
        self.func_rayleigh = RegularGridInterpolator(
            (self.WL, self.AIRMASS), data_rayleigh
        )
        self.func_O2abs = RegularGridInterpolator((self.WL, self.AIRMASS), data_O2abs)
        self.func_PWVabs = RegularGridInterpolator(
            (self.WL, self.AIRMASS, self.PWV), data_PWVabs
        )
        self.func_OZabs = RegularGridInterpolator(
            (self.WL, self.AIRMASS, self.OZ), data_OZabs
        )

    def __str__(self):
        return self.Name

    #
    # functions to access to interpolated transparency functions on training dataset
    #
    def GetWL(self):
        """Return wavelength array used by the grid"""
        return self.WL

    def GetRayleighTransparencyArray(self, wl, am):
        """Return Rayleigh transmission for the corresponding wavelength array
        at the airmass"""
        pts = [(the_wl, am) for the_wl in wl]
        pts = np.array(pts)
        return self.func_rayleigh(pts)

    def GetO2absTransparencyArray(self, wl, am):
        """Return O2 transmission for the corresponding wavelength array
        at the airmass"""
        pts = [(the_wl, am) for the_wl in wl]
        pts = np.array(pts)
        return self.func_O2abs(pts)

    def GetPWVabsTransparencyArray(self, wl, am, pwv):
        """Return PWV transmission for the corresponding wavelength array
        at the airmass"""
        pts = [(the_wl, am, pwv) for the_wl in wl]
        pts = np.array(pts)
        return self.func_PWVabs(pts)

    def GetOZabsTransparencyArray(self, wl, am, oz):
        """Return Ozone transmission for the corresponding wavelength array
        at the airmass"""
        pts = [(the_wl, am, oz) for the_wl in wl]
        pts = np.array(pts)
        return self.func_OZabs(pts)

    def GetGriddedTransparencies(
        self,
        wl,
        am,
        pwv,
        oz,
        flagRayleigh=True,
        flagO2abs=True,
        flagPWVabs=True,
        flagOZabs=True,
    ):
        """
        Emulation of libradtran simulated transparencies. Decomposition of the
        total transmission in different processes:
        - Rayleigh scattering
        - O2 absorption
        - PWV absorption
        - Ozone absorption
        Those 4 processes transmission are derived from an interpolated grid
        of the atmospheric parameters arguments.

        Note the aerosols are not included as their modelisation is outside

        :param wl: wavelength array or list
        :type wl: float

        :param am: the airmass
        :type am: float from 1.0 to 2.5

        :param pwv: the precipitable water vapor
        :type pwv: float in unit of mm

        :param oz: the ozone column depth
        :type oz: float in Dobson unit

        :param flagRayleigh: flags to activate Rayleigh scattering process,
        default True
        :type flagRayleigh: bool, optional

        :param flagO2abs: flags to activate Oxygen absorption process,
        default True
        :type flagO2abs: bool, optional

        :param flagPWVabs: flags to activate Precipitable water vapor
        absorption process, default True
        :type flagPWVabs: bool, optional

        :param flagOZabs: flags to activate Ozone absorption process,
        default True
        :type flagOZabs: bool, optional

        :return: array of atmospheric transmissions (save size as `wl`)
        :rtype: floats

        """

        if flagRayleigh:
            transm = self.GetRayleighTransparencyArray(wl, am)
        else:
            transm = np.ones(len(wl))

        if flagO2abs:
            transm *= self.GetO2absTransparencyArray(wl, am)

        if flagPWVabs:
            transm *= self.GetPWVabsTransparencyArray(wl, am, pwv)

        if flagOZabs:
            transm *= self.GetOZabsTransparencyArray(wl, am, oz)

        return transm

    def GetAerosolsTransparencies(self, wl, am, tau=0.0, beta=1.4):
        """
        Compute transmission due to aerosols.
        A simple model of aerosols is assumed based on one component with one
        optical depth and one Angstrom exponent.

        :param wl: wavelength array or list
        :type wl: float

        :param am: the airmass
        :type am: float from 1.0 to 2.5

        :param tau: the vertical aerosol depth of each component at lambda0
        wavelength, default set to 0.0 for no aerosol component
        :type tau: float

        :param beta: the angstrom exponent. Must be positive in the range
        0.0 to 3.0
        :type beta: float

        :return: 1D array of atmospheric transmission (save size as wl)
        :rtype: array of floats

        """

        wl = np.array(wl)

        exponent = (tau / self.tau0) * np.exp(-beta * np.log(wl / self.lambda0)) * am

        transm = np.exp(-exponent)
        return transm

    def GetAllTransparencies(
        self,
        wl,
        am,
        pwv,
        oz,
        tau=0.0,
        beta=1.4,
        flagRayleigh=True,
        flagO2abs=True,
        flagPWVabs=True,
        flagOZabs=True,
        flagAerosols=True,
    ):
        """
        Combine interpolated libradtran transmission with analytical expression
        for the aerosols

        :param wl: wavelength array or list
        :type wl: float

        :param am: the airmass
        :type am: float from 1.0 to 2.5

        :param pwv: the precipitable water vapor
        :type pwv: float in unit of mm

        :param oz: the ozone column depth
        :type oz: float in Dobson unit

        :param tau: the vertical aerosol depth of each component at lambda0
        wavelength, default set to 0.0 for no aerosol component
        :type tau: float

        :param beta: the angstrom exponent. Must be positive in the range 0., 3.
        :type beta: float

        :param flagRayleigh: flags to activate Rayleigh scattering process,
        default True
        :type flagRayleigh: bool, optional

        :param flagO2abs: flags to activate Oxygen absorption process,
        default True
        :type flagO2abs: bool, optional

        :param flagPWVabs: flags to activate Precipitable water vapor
        absorption process, default True
        :type flagPWVabs: bool, optional

        :param flagOZabs: flags to activate Ozone absorption process,
        default True
        :type flagOZabs: bool, optional

        :param flagAerosols: flags to activate Aerosol scattering,
        default True
        :type flagAerosols: bool, optional

        :return: 1D array of atmospheric transmission (save size as wl)
        :rtype: array of floats


        """

        transm = self.GetGriddedTransparencies(
            wl,
            am,
            pwv,
            oz,
            flagRayleigh=flagRayleigh,
            flagO2abs=flagO2abs,
            flagPWVabs=flagPWVabs,
            flagOZabs=flagOZabs,
        )

        if flagAerosols:
            transmaer = self.GetAerosolsTransparencies(wl, am, tau, beta)
            transm *= transmaer

        return transm


class ObsAtmoPressure(ObsAtmoGrid):
    """
    Emulate Atmospheric Transparency above LSST from a data grids
    extracted from libradtran and analytical functions for aerosols.
    There are 3 grids:
    - 2D grid Rayleigh transmission vs (wavelength,airmass)
    - 2D grid O2 absorption vs  (wavelength,airmass)
    - 3D grid for PWV absorption vs (wavelength,airmass,PWV)
    - 3D grid for Ozone absorption vs (wavelength,airmass,Ozone)
    - Aerosol transmission for only one component.

    It extend the functionalities of ObsAtmoGrid Emulator by providing a
    correction when the ground pressure differs from the standard pressure
    expected at that site.

    :param obs_label:  The observation site name. Five sites have been
    preselected.
    :type obs_label: string among 'LSST','CTIO','OHP','PDM','OMK','OSL'

    :param pressure: the local pressure. Set it to zero (sefault value)
    if the standard pressure expected must be used.
    :type pressure: float in unit of hPa or mbar. Optional.

    :raise: the `obs_label` is not in the list of predefined observatory sites

    :return: the emulator obsject
    :rtype: object of class `ObsAtmoPressure`

    Examples
    --------
    >>> emulator = ObsAtmoPressure('LSST')

    """

    def __init__(self, obs_str="LSST", pressure=0):
        """
        Initialize the class for data point files from which the 2D and 3D grids
        are created.
        Interpolation are calculated from the scipy RegularGridInterpolator()
        function
        Both types of data : training data for normal interpolation use and
        the test data used to check accuracy of the interpolation of data.

        :param obs_label:  The observation site name. Five sites have been
        preselected.
        :type obs_label: string among 'LSST','CTIO','OHP','PDM','OMK','OSL'

        :param pressure: the local pressure. Set it to zero (sefault value)
        if the standard pressure
        expected must be used.
        :type pressure: float in unit of hPa or mbar. Optional.

        :raises: the `obs_label` is not in the list of predefined
        observatory sites

        :return: the emulator obsject
        :rtype: object of class `ObsAtmoPressure`

        """
        ObsAtmoGrid.__init__(self, obs_str=obs_str)

        self.pressure = pressure
        # self.refpressure = Dict_Of_sitesPressures[obs_str]
        self.refpressure = Dict_Of_sitesPressures[self.OBS_tag]

        if self.refpressure <= 0:
            raise ValueError(
                f"Invalid reference pressure for {self.OBS_tag}:" f" {self.refpressure}"
            )

        self.pressureratio = self.pressure / self.refpressure
        if pressure == 0.0:
            self.pressureratio = 1
            self.pressure = self.refpressure

        self.Name = "Atmospheric emulator ObsAtmoPressure for "
        self.Name = (
            self.Name + f"observation site {self.OBS_tag} P = {self.pressure} hPa"
        )

    def GetRayleighTransparencyArray(self, wl, am):
        """
        Correct the Rayleigh scattering implemented in the base class
        `ObsAtmoGrid`, because the pressure is different from the standard one.

        :param wl: wavelength array or list
        :type wl: float

        :param am: the airmass
        :type am: float from 1.0 to 2.5

        :return: array of air transparencies due to Rayleigh scattering
        at the corresponding pressure
        :rtype: array of floats

        """
        return np.power(
            super().GetRayleighTransparencyArray(wl, am), self.pressureratio
        )

    def GetO2absTransparencyArray(self, wl, am, satpower=1.16306918):
        """
        Correction of O2 absorption profile by the P/Pref with a power
        estimated from libradtran simulations, where P is the true pressure
        and Pref is the reference pressure for the site.

        :param wl: wavelength array or list
        :type wl: float

        :param am: the airmass
        :type am: float from 1.0 to 2.5

        :return: array of air transparencies due to Oxygen absorption
        at the corresponding pressure
        :rtype: array of floats

        """
        return np.power(
            super().GetO2absTransparencyArray(wl, am),
            np.power(self.pressureratio, satpower),
        )


class ObsAtmo(ObsAtmoPressure):
    """
    Emulate Atmospheric Transparency above different sites.
    The preselected sites are LSST,CTIO, Mauna Kea, Observatoire de
    Haute Provence, Pic du Midi or Sea Level.
    Each site corresponds has an corresponding pressure.
    If the pressure does not correspond
    to the standard one, the pressure can be renormalized.

    ObsAtmo is the user end-point which call the official implementation of the
    emulator.

    By now ``ObsAtmo`` refer to the ``ObsAtmoPressure`` which itself rely on the
    ``ObsAtmoGrid``.

    :param obs_label:  The observation site name. Five sites have been
    preselected.
    :type obs_label: string among 'LSST','CTIO','OHP','PDM','OMK','OSL'

    :param pressure: the local pressure. Set it to zero (sefault value)
    if the standard pressure
    expected must be used.
    :type pressure: float in unit of hPa or mbar. Optional.

    :raises: the `obs_label` is not in the list of predefined
    observatory sites

    :return: the emulator obsject
    :rtype: object of class `ObsAtmo`


    Usage
    -----

    .. doctest::

        >>> emul = ObsAtmo('LSST', 743.0)
        >>> wl = [400., 800., 900.]
        >>> am=1.2
        >>> pwv=4.0
        >>> oz=300.
        >>> transm = emul.GetAllTransparencies(wl,am,pwv,oz)
        >>> print(wl)
        [400.0, 800.0, 900.0]
        >>> print(transm)
        [0.72485491 0.97330618 0.85675228]
    """

    def __init__(self, obs_str="LSST", pressure=0):
        """
        Initialize the ``ObsAtmo``

        :param obs_label:  The observation site name. Five sites have been
        preselected.
        :type obs_label: string among 'LSST','CTIO','OHP','PDM','OMK','OSL'

        :param pressure: the local pressure. Set it to zero (sefault value)
        if the standard pressure
        expected must be used.
        :type pressure: float in unit of hPa or mbar. Optional.

        :raise: the `obs_label` is not in the list of predefined
        observatory sites

        :return: the emulator obsject
        :rtype: object of class `ObsAtmo`


        Examples
        --------

        .. doctest::

            >>> e = ObsAtmo(obs_str="LSST", pressure=0)
            >>> print(e)
            Atmospheric emulator ObsAtmo for observation site LSST

        """
        ObsAtmoPressure.__init__(self, obs_str=obs_str, pressure=pressure)
        # self.Name = f"Atmospheric emulator ObsAtmo for observation site {obs_str}"
        self.Name = "Atmospheric emulator ObsAtmo for observation site "
        self.Name = self.Name + f"{self.OBS_tag}"

    def plot_transmission(
        self,
        am=1.0,
        pwv=4.0,
        oz=400.0,
        tau=0.1,
        beta=1.4,
        xscale="linear",
        yscale="linear",
        savepath=None,
    ):
        """Plot ObsAtmo transmission

        :param am: the airmass, default set to 1.0
        :type am: float from 1.0 to 2.5, optional

        :param pwv: the precipitable water vapor, default set to 4.0 mm
        :type pwv: float in unit of mm, optional

        :param oz: the ozone column depth, default set to 400 DU.
        :type oz: float in Dobson unit, optional

        :param tau: the vertical aerosol depth of each component at lambda0
        wavelength,default set to 0.1
        :type tau: float,optional

        :param beta: the angstrom exponent. Must be positive in the
        range 0-3.0. Default set to +1.2
        :type beta: float, optional

        Examples
        --------
        >>> e = ObsAtmo(obs_str = "LSST", pressure = 0)
        >>> e.plot_transmission()

        """
        wls = self.WL
        transm = self.GetAllTransparencies(wls, am, pwv, oz, tau, beta)

        textstr = "\n".join(
            (
                f"airmass = {am:.1f}",
                f"pwv = {pwv:.1f} mm",
                f"ozone = {oz:.0f} DU",
                f"$\\tau$ = {tau:.3f}",
                f"$\\beta$ = {beta:.1f}",
            )
        )
        props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)

        fig, ax = plt.subplots()
        ax.plot(wls, transm, "b-")
        ax.grid()
        ax.set_yscale(yscale)
        ax.set_xscale(xscale)
        ax.set_title(
            f"Atmospheric transmission at {self.OBS} with "
            f" P={self.pressure:.1f} hPa"
        )
        ax.set_xlabel(r"$\lambda$ [nm]")
        ax.set_ylabel("Transmission")
        # place a text box in upper left in axes coords
        ax.text(
            0.75,
            0.05,
            textstr,
            transform=ax.transAxes,
            fontsize=14,
            verticalalignment="bottom",
            bbox=props,
        )

        if savepath:
            fig.savefig(savepath, dpi=200)
        else:
            fig.show()


def usage():
    print("*******************************************************************")
    print(sys.argv[0], " -s<observation site-string> -p pressure")
    print("Observation sites are : ")
    print(" ".join(Dict_Of_sitesAltitudes.keys()))
    print("if pressure is not given", "the standard pressure for the site is used")

    print("\t actually provided : ")
    print("\t \t Number of arguments:", len(sys.argv), "arguments.")
    print("\t \t Argument List:", str(sys.argv))


def run(obs_str, pressure):
    print("==================================================================")
    print(
        f"Atmospheric Emulator for {obs_str} observatory and "
        f"pressure = {pressure:.2f} hPa"
    )
    print("==================================================================")

    emul = ObsAtmo(obs_str=obs_str, pressure=pressure)
    wl = [400.0, 800.0, 900.0]
    am = 1.2
    pwv = 4.0
    oz = 300.0
    transm = emul.GetAllTransparencies(wl, am, pwv, oz)
    print("wavelengths (nm) \t = ", wl)
    print("transmissions    \t = ", transm)


def is_float(element: Any) -> bool:
    # If you expect None to be passed:
    if element is None:
        return False
    try:
        float(element)
        return True
    except ValueError:
        return False


if __name__ == "__main__":

    import doctest
    import getopt

    doctest.testmod()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hs:p:", ["s=", "p="])
    except getopt.GetoptError:
        raise ValueError(
            "Exception bad getopt with :: "
            + sys.argv[0]
            + " -s<observation-site-string>"
        )

    print("opts = ", opts)
    print("args = ", args)

    obs_str = ""
    pressure_str = ""

    for opt, arg in opts:
        if opt == "-h":
            usage()
            sys.exit(0)
        elif opt in ("-s", "--site"):
            obs_str = arg
        elif opt in ("-p", "--pressure"):
            pressure_str = arg

    if is_float(pressure_str):
        pressure = float(pressure_str)
    elif pressure_str == "":
        pressure = 0
    else:
        print(f"Pressure argument {pressure_str} is not a float")
        sys.exit()

    if obs_str in Dict_Of_sitesAltitudes.keys():
        run(obs_str=obs_str, pressure=pressure)
    else:
        raise ValueError(
            f"Observatory {obs_str} not in preselected observation site.\n"
            f"This site {obs_str} must be added in libradtranpy "
            f"preselected sites."
        )
