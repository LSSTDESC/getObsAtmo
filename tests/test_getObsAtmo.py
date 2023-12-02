# python -m unittest
# python -m unittest -v

import unittest
from getObsAtmo import Dict_Of_sitesAltitudes,Dict_Of_sitesAliases,Dict_Of_sitesPressures,file_data_dict
from getObsAtmo import _getPackageDir,sanitizeString,validateObsName,is_obssite, ObsAtmo
import numpy as np 
import os



class GetObsAtmoTestCase(unittest.TestCase):
    """A test case for the getObsAtmo package."""

    def test_packagedir(self):
        """Test if the installation dir and data dir exists
        """
        installationdir = _getPackageDir()
        # construct the path of input data files
        datapath = os.path.join(installationdir, '../obsatmo_data')
        self.assertTrue(os.path.exists(datapath))
        self.assertTrue(os.path.isdir(datapath))    

    def test_inputdata(self):
        """Test if input data exists
        """
        installationdir = _getPackageDir()
        datapath = os.path.join(installationdir, '../obsatmo_data')
        for obs_str in Dict_Of_sitesAltitudes.keys():
            for file_key,filenamepart in file_data_dict.items():
                filename = obs_str + "_" + filenamepart
                fullfilename = os.path.join(datapath,filename)
                self.assertTrue(os.path.exists(fullfilename))
                self.assertTrue(os.path.isfile(fullfilename))

    def test_sanitizedstring(self):
        """Test sanitized string
        """
        self.assertEqual(sanitizeString('L s s T'),'LSST')

    def test_validateobsname(self):
        """Test validate observatory name
        """
        self.assertEqual(validateObsName('r u bin'),'LSST')
        self.assertEqual(validateObsName('rubin observatory'),'LSST')
        self.assertEqual(validateObsName('A u x Tel'),'LSST')
        self.assertEqual(validateObsName('cerro tololo'),'CTIO')
        self.assertEqual(validateObsName('pic du midi'),'PDM')
        self.assertEqual(validateObsName('Observatoire du Pic du Midi'),'PDM')
        self.assertEqual(validateObsName('mauna kea'),'OMK')
        self.assertEqual(validateObsName('mauna kea observatory'),'OMK')
        self.assertEqual(validateObsName('sea level'),'OSL')
        self.assertEqual(validateObsName('sea level observatory'),'OSL')
        self.assertIsNone(validateObsName('AaBb'))


    def test_is_obssite(self):
        """Test if observation site exist
        """
        self.assertTrue(is_obssite("LSST"))
        self.assertTrue(is_obssite("a u x tel"))
        self.assertTrue(is_obssite("AUXTEL"))
        self.assertTrue(is_obssite("CTIO"))
        self.assertTrue(is_obssite("OMK"))
        self.assertFalse(is_obssite("AMissingObs"))
        self.assertTrue(is_obssite("L s s t"))

    def test_ObsAtmo(self):
        """Test ObsAtmo emulator
        """
        e = ObsAtmo('LSST')
        array = e.GetAllTransparencies([300.,400.,600.,800.,1000.],1.,0.,0.)
        self.assertTrue(np.allclose(array,[0.41483073, 0.76797745, 0.95170851, 0.98476396, 0.9937721 ]),msg="ObsAtmo.GetAllTransparencies error")

        e = ObsAtmo('AUXTEL')
        array = e.GetAllTransparencies([300.,400.,600.,800.,1000.],1.,0.,0.)
        self.assertTrue(np.allclose(array,[0.41483073, 0.76797745, 0.95170851, 0.98476396, 0.9937721 ]),msg="ObsAtmo.GetAllTransparencies error")


if __name__ == "__main__":
    unittest.main()
