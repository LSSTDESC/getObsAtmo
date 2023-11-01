import unittest
from getObsAtmo import is_obssite, ObsAtmo
import numpy as np 



class GetObsAtmoTestCase(unittest.TestCase):
    """A test case for the getCalspec package."""

    def test_is_obssite(self):
        self.assertTrue(is_obssite("LSST"))
        self.assertTrue(is_obssite("CTIO"))
        self.assertTrue(is_obssite("OMK"))
        self.assertFalse(is_obssite("AMissingObs"))
        self.assertTrue(is_obssite("L s s t"))

    def test_ObsAtmo(self):
        e = ObsAtmo('LSST')
        array = e.GetAllTransparencies([300.,400.,600.,800.,1000.],1.,0.,0.)
        self.assertTrue(np.allclose(array,[0.41483073, 0.76797745, 0.95170851, 0.98476396, 0.9937721 ]),msg="ObsAtmo.GetAllTransparencies error")

        


if __name__ == "__main__":
    unittest.main()
