import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pickle
import matplotlib.pyplot as plt
import pandas as pd
import os
import warnings
import sys,getopt



__all__ = ['get_obssite_keys',
           'is_obssite',
           'ObsAtmo',
           ]


# preselected sites 
Dict_Of_sitesAltitudes = {'LSST':2.663,
                          'CTIO':2.207,
                          'OHP':0.65,
                          'PDM':2.8905,
                          'OMK':4.205,
                          'OSL':0.000,
                           }
# pressure calculated by libradtran
Dict_Of_sitesPressures = {'LSST':731.50433,
                          'CTIO':774.6052,
                          'OHP':937.22595,
                          'PDM':710.90637,
                          'OMK':600.17224,
                          'OSL':1013.000,
                        }

file_data_dict = {
    "info" :"atmospherictransparencygrid_params.pickle",
    "data_rayleigh" : "atmospherictransparencygrid_rayleigh.npy",
    "data_o2abs" : "atmospherictransparencygrid_O2abs.npy",
    "data_pwvabs": "atmospherictransparencygrid_PWVabs.npy",
    "data_ozabs" : "atmospherictransparencygrid_OZabs.npy",
}

def getObsSiteDataFrame():  
    df = pd.DataFrame(columns=['altitude','pressure'], index=list(Dict_Of_sitesAltitudes.keys()))
    for key in Dict_Of_sitesAltitudes.keys():
        df.loc[key] = pd.Series({'altitude':Dict_Of_sitesAltitudes[key],'pressure':Dict_Of_sitesPressures[key]})
    return df


def _getPackageDir():
    """This method must live in the top level of this package, so if this
    moves to a utils file then the returned path will need to account for that.
    """
    dirname = os.path.dirname(__file__)
    return dirname


def sanitizeString(label):
    """This method sanitizes the star label."""
    return label.upper().replace(' ','')


def sanitizeDataFrame(df):
    """This method sanitizes the star label."""
    tmp_df = df.str.upper()
    tmp_df = tmp_df.str.replace(' ', '')
    return tmp_df


def get_obssite_keys(obs_label):
    """Return the DataFrame keys if a star name corresponds to a Calspec entry
    in the tables.

    Parameters
    ----------
    obs_label: str
        The observation site name.

    Returns
    -------
    keys: array_like
        The DataFrame keys corresponding to the star name.

    Examples
    --------
    >>> get_obssite_keys("lsst")   #doctest: +ELLIPSIS
    0      False
    1      False
    ...
    """
    label = sanitizeString(obs_label)
    df = getObsSiteDataFrame()
    name_index = [name.upper()  for name in df.index]
    if len(name_index) > 0:
        keys = pd.Series([False] * len(df))
        for name in name_index:
            tmp_df = sanitizeDataFrame(df[name])
            keys = keys | (tmp_df == label)
        return keys
    else:
        raise KeyError("No observation site label {obs_label} in config dictionaries")


def is_obssite(obs_label):
    """
    Test if an observation site label is among the list atmospheric grids.

    Parameters
    ----------
    obs_label: str
        The observation site name.

    Returns
    -------
    is_obssite: bool
        True if the observation site is star has a table.

    Examples
    --------
    >>> is_obssite("LSST")
    True
    >>> is_obssite("OHP")
    True
    """
    return np.any(get_obssite_keys(sanitizeString(obs_label)))

 
class ObsAtmoGrid:
    """"
    Emulate Atmospheric Transparency above LSST from a data grids
    extracted from libradtran and analytical functions for aerosols.
    There are 3 grids:
    - 2D grid Rayleigh transmission vs (wavelength,airmass)
    - 2D grid O2 absorption vs  (wavelength,airmass)
    - 3D grid for PWV absorption vs (wavelength,airmass,PWV)
    - 3D grid for Ozone absorption vs (wavelength,airmass,Ozone)
    - Aerosol transmission for any number of components
    """
    def __init__(self,obs_str = "LSST") :
        """
        Initialize the class for data point files from which the 2D and 3D grids are created.
        Interpolation are calculated from the scipy RegularGridInterpolator() function
        Both types of data : trainging data for normal interpolaton use and the test data used
        to check accuracy of the interpolation of data.

        Parameters
        ----------
            obs_str : str
                pre-defined observation site tag corresponding to data files in data path
            
        Returns
        --------
            the emulator


        """
        OBS_tag = ""
        if obs_str in Dict_Of_sitesAltitudes.keys():
            OBS_tag = obs_str
            print(f"Observatory {obs_str} found in preselected observation sites")
        else:
            print(f"Observatory {obs_str} not in preselected observation sites")
            print(f"This site {obs_str} must be added in libradtranpy preselected sites")
            print(f"and generate corresponding scattering and absorption profiles")
            sys.exit()

        
        self.Name = f"Atmospheric emulator ObsAtmoGrid for observation site {OBS_tag}"

        # construct the path of input data files
        self.path = os.path.join(_getPackageDir(),'../obsatmo_data')
        self.fn_info = OBS_tag + "_" + file_data_dict["info"]
        self.fn_rayleigh = OBS_tag + "_" + file_data_dict["data_rayleigh"]
        self.fn_O2abs = OBS_tag + "_" + file_data_dict["data_o2abs"]
        self.fn_PWVabs = OBS_tag + "_" + file_data_dict["data_pwvabs"]
        self.fn_OZabs = OBS_tag + "_" + file_data_dict["data_ozabs"]
       
        self.info_params = None
        self.data_rayleigh = None
        self.data_O2abs = None
        self.data_PWVabs = None
        self.data_OZabs = None
     
        
        # load all data files (training and test)
        self.loadtables()
        
        # setup training dataset (those used for interpolation)
        self.WLMIN = self.info_params["WLMIN"]
        self.WLMAX = self.info_params["WLMAX"]
        self.WLBIN = self.info_params["WLBIN"]
        self.NWLBIN = self.info_params['NWLBIN']
        self.WL = self.info_params['WL']
        self.OBS = self.info_params['OBS']
        
        self.AIRMASSMIN = self.info_params['AIRMASSMIN']
        self.AIRMASSMAX = self.info_params['AIRMASSMAX']
        self.NAIRMASS = self.info_params['NAIRMASS']
        self.DAIRMASS = self.info_params['DAIRMASS']
        self.AIRMASS = self.info_params['AIRMASS']
        
        self.PWVMIN = self.info_params['PWVMIN']
        self.PWVMAX = self.info_params['PWVMAX'] 
        self.NPWV = self.info_params['NPWV']
        self.DPWV = self.info_params['DPWV'] 
        self.PWV = self.info_params['PWV']
        
        
        self.OZMIN =  self.info_params['OZMIN']
        self.OZMAX = self.info_params['OZMAX']
        self.NOZ = self.info_params['NOZ']
        self.DOZ =  self.info_params['DOZ'] 
        self.OZ = self.info_params['OZ']
        

        # constant parameters defined for aerosol formula
        self.lambda0 = 550.
        self.tau0 = 1.

        # interpolation is done over training dataset
        self.func_rayleigh = RegularGridInterpolator((self.WL,self.AIRMASS),self.data_rayleigh)
        self.func_O2abs = RegularGridInterpolator((self.WL,self.AIRMASS),self.data_O2abs)
        self.func_PWVabs = RegularGridInterpolator((self.WL,self.AIRMASS,self.PWV),self.data_PWVabs)
        self.func_OZabs = RegularGridInterpolator((self.WL,self.AIRMASS,self.OZ),self.data_OZabs)

        
        
    def loadtables(self):
        """
        Load files into grid arrays
        """
        
        filename=os.path.join(self.path,self.fn_info)     
        with open(filename, 'rb') as f:
            self.info_params = pickle.load(f)
            
       
        filename=os.path.join(self.path,self.fn_rayleigh)
        with open(filename, 'rb') as f:
            self.data_rayleigh = np.load(f)
            
            
        filename=os.path.join(self.path,self.fn_O2abs)
        with open(filename, 'rb') as f:
            self.data_O2abs = np.load(f)
            
      
        filename=os.path.join(self.path,self.fn_PWVabs)
        with open(filename, 'rb') as f:
            self.data_PWVabs = np.load(f)
            
        
            
        filename=os.path.join(self.path,self.fn_OZabs)
        with open(filename, 'rb') as f:
            self.data_OZabs = np.load(f)
            

    def __str__(self):
        return self.Name
     
    #        
    # functions to access to interpolated transparency functions on training dataset
    #        
    def GetWL(self):
        return self.WL
    
    def GetRayleighTransparencyArray(self,wl,am):
        pts = [ (the_wl,am) for the_wl in wl ]
        pts = np.array(pts)
        return self.func_rayleigh(pts)
    
    
    def GetO2absTransparencyArray(self,wl,am):
        pts = [ (the_wl,am) for the_wl in wl ]
        pts = np.array(pts)
        return self.func_O2abs(pts)
    
    
    def GetPWVabsTransparencyArray(self,wl,am,pwv):
        pts = [ (the_wl,am,pwv) for the_wl in wl ]
        pts = np.array(pts)
        return self.func_PWVabs(pts)
    
    
    def GetOZabsTransparencyArray(self,wl,am,oz):
        pts = [ (the_wl,am,oz) for the_wl in wl ]
        pts = np.array(pts)
        return self.func_OZabs(pts)
    
        
    def GetGriddedTransparencies(self,wl,am,pwv,oz,flagRayleigh=True,flagO2abs=True,flagPWVabs=True,flagOZabs=True):
        """
        Emulation of libradtran simulated transparencies. Decomposition of the
        total transmission in different processes:
        - Rayleigh scattering
        - O2 absorption
        - PWV absorption
        - Ozone absorption
        
        Parameters
        ----------
          - wl : wavelength array or list
          - am :the airmass,
          - pwv : the precipitable water vapor (mm)
          - oz : the ozone column depth in Dobson unit
          - flags to activate or not the individual interaction processes
        
        Returns
        ------
           - 1D array of atmospheric transmission (save size as wl)
        
        """
        

        if flagRayleigh:
            transm = self.GetRayleighTransparencyArray(wl,am)
        else:
            transm = np.ones(len(wl))
            
        if flagO2abs:
            transm *= self.GetO2absTransparencyArray(wl,am)
            
        if flagPWVabs:
            transm *= self.GetPWVabsTransparencyArray(wl,am,pwv)
            
        if flagOZabs:
            transm *= self.GetOZabsTransparencyArray(wl,am,oz)
            
        return transm
            
    def GetAerosolsTransparencies(self,wl,am,tau=0,beta=-1):
        """
        Compute transmission due to aerosols:
        
        Parameters
        ----------
            - wl : wavelength array
            - am : the airmass
            - tau : the vertical aerosol depth of each component at lambda0 vavelength
            - beta : the angstrom exponent. Must be negativ.
        
        Returns
        -------
            - 1D array of atmospheric transmission (save size as wl)
        
        """
          
        wl = np.array(wl)
                
        exponent = (tau/self.tau0)*np.exp(beta*np.log(wl/self.lambda0))*am
        transm = np.exp(-exponent)
            
        return transm
        
        
    def GetAllTransparencies(self,wl,am,pwv,oz, tau=0, beta=-1, flagRayleigh=True,flagO2abs=True,flagPWVabs=True,flagOZabs=True,flagAerosols=False):
        """
        Combine interpolated libradtran transmission with analytical expression for the
        aerosols
        
        Parameters
        ----------

            - wl : wavelength array or list
            - am :the airmass,
            - pwv : the precipitable water vapor (mm)
            - oz : the ozone column depth in Dobson unit
            - tau & beta : parameters for aerosols
            - flags to activate or not the individual interaction processes
        
        Returns
        -------
            - 1D array of atmospheric transmission (save size as wl)
        
        """
        
        transm = self.GetGriddedTransparencies(wl,am,pwv,oz,flagRayleigh=flagRayleigh,flagO2abs=flagO2abs,flagPWVabs=flagPWVabs,flagOZabs=flagOZabs)
        
        if flagAerosols:
            transmaer = self.GetAerosolsTransparencies(wl,am,tau,beta)
            transm *=transmaer
           
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
    - Aerosol transmission for any number of components

    It uses the SimpleAtm Emulator. This particular class interpolate transparency
    with local pressures.

    """
    def __init__(self,obs_str = "LSST", pressure = 0 ) : 
        ObsAtmoGrid.__init__(self,obs_str = obs_str)
        """
        Initialize the class for data point files from which the 2D and 3D grids are created.
        Interpolation are calculated from the scipy RegularGridInterpolator() function
        Both types of data : trainging data for normal interpolaton use and the test data used
        to check accuracy of the interpolation of data.

        Parameters
        ----------

            obs_str : pre-defined observation site tag corresponding to data files in data path
            pressure : pressure for which one want the transmission in mbar or hPa
         
        Returns
        -------
            the emulator object ObsAtmoPressure

        """
    
        self.pressure = pressure
        self.refpressure = Dict_Of_sitesPressures[obs_str]
        self.pressureratio = self.pressure/self.refpressure
        if pressure == 0.0:
            self.pressureratio = 1


        self.Name = f"Atmospheric emulator ObsAtmoPressure for observation site {obs_str}"

    def GetRayleighTransparencyArray(self,wl,am):
        """
        Scaling of optical depth by the term P/Pref, where P is the true pressure
        and Pref is the reference pressure for the site.
        """
        return np.power(super().GetRayleighTransparencyArray(wl,am),self.pressureratio)
    

    def GetO2absTransparencyArray(self,wl,am,satpower=1.16306918):
        """
        Correction of O2 absorption profile by the P/Pref with a power estimated
        from libradtran simulations, where P is the true pressure
        and Pref is the reference pressure for the site.

        Comparing LSST site with pressure at Mauna Kea and Sea Level show the satpower
        = 1.1 is appropriate.
        """
        return np.power(super().GetO2absTransparencyArray(wl,am),
                        np.power(self.pressureratio,satpower))


class ObsAtmo(ObsAtmoPressure):
    
    """
    Emulate Atmospheric Transparency above different sites.
    The preselected sites are LSST,CTIO, Mauna Kea, Observatoire de Haute Provence,
    Pic du Midi or Sea Level.
    Each site corresponds has an corresponding pressure. If the pressure does not correspond
    to the standard one, the pressure can be renormalized.

    ObsAtmo is the user end-point which call the official implementation of the emulator.
    
    By now ``ObsAtmo`` refer to the ``ObsAtmoPressure`` which itself rely on the ``ObsAtmoGrid``.
    """
    def __init__(self,obs_str = "LSST", pressure = 0 ) : 
        ObsAtmoPressure.__init__(self,obs_str = obs_str, pressure = pressure )
        """
        Initialize the ``ObsAtmo`` 

        Parameters 
        ----------
          obs_str : str 
               pre-defined observation site tag corresponding to data files in data path
          pressure : float
               pressure for which one want the transmission in mbar or hPa

        Returns
        -------
            the emulator object ObsAtmoPressure
        
        Examples
        --------
        >>> e = ObsAtmo(obs_str = "LSST", pressure = 0)
        >>> print(c)   #doctest: +ELLIPSIS


        """
        self.Name = f"Atmospheric emulator ObsAtmo for observation site {obs_str}"



def usage():
    print("*******************************************************************")
    print(sys.argv[0],' -s<observation site-string> -p pressure')
    print("Observation sites are : ")
    print(' '.join(Dict_Of_sitesAltitudes.keys()))
    print('if pressure is not given, the standard pressure for the site is used')


    print('\t actually provided : ')
    print('\t \t Number of arguments:', len(sys.argv), 'arguments.')
    print('\t \t Argument List:', str(sys.argv))


def run(obs_str, pressure):
    print("===========================================================================")
    print(f"Atmospheric Emulator for {obs_str} observatory and pressure = {pressure:.2f} hPa") 
    print("===========================================================================")
    
    
    emul = ObsAtmo(obs_str = obs_str, pressure = pressure)
    wl = [400.,800.,900.]
    am=1.2
    pwv =4.0
    oz=300.
    transm = emul.GetAllTransparencies(wl,am,pwv,oz)
    print("wavelengths (nm) \t = ",wl)
    print("transmissions    \t = ",transm)

def is_float(element: any) -> bool:
    #If you expect None to be passed:
    if element is None: 
        return False
    try:
        float(element)
        return True
    except ValueError:
        return False

if __name__ == "__main__":

    import doctest

    doctest.testmod()

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hs:p:",["s=","p="])
    except getopt.GetoptError:
        print(' Exception bad getopt with :: '+sys.argv[0]+ ' -s<observation-site-string>')
        sys.exit(2)

    print('opts = ',opts)
    print('args = ',args)

    obs_str = ""
    pressure_str =""

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-s", "--site"):
            obs_str = arg
        elif opt in ("-p", "--pressure"):
            pressure_str = arg
       
    if is_float(pressure_str):
        pressure = float(pressure_str)
    elif pressure_str =="":
        pressure = 0
    else:
        print(f"Pressure argument {pressure_str} is not a float")
        sys.exit()

        
    if obs_str in Dict_Of_sitesAltitudes.keys():
        run(obs_str=obs_str,pressure = pressure)
    else:
        print(f"Observatory {obs_str} not in preselected observation site")
        print(f"This site {obs_str} must be added in libradtranpy preselected sites")
        sys.exit()


