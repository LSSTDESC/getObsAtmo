import unittest
import numpy as np
import os
from getObsAtmo import (
    Dict_Of_sitesAltitudes,
    file_data_dict,
    _getPackageDir,
    sanitizeString,
    validateObsName,
    is_obssite,
    ObsAtmo,
    ObsAtmoPressure,
    ObsAtmoGrid,
    convert_dict_to_json,
    convert_json_to_dict,
)


class GetObsAtmoTestCase(unittest.TestCase):
    """A test case for the getObsAtmo package."""

    def test_packagedir(self):
        """Test if the installation dir and data dir exists"""
        installationdir = _getPackageDir()
        datapath = os.path.join(installationdir, "obsatmo_data")

        self.assertTrue(os.path.exists(datapath))
        self.assertTrue(os.path.isdir(datapath))

    def test_inputdata(self):
        """Test if input data exists"""
        installationdir = _getPackageDir()
        datapath = os.path.join(installationdir, "obsatmo_data")
        for obs_str in Dict_Of_sitesAltitudes.keys():
            for file_key, filenamepart in file_data_dict.items():
                if file_key == obs_str:
                    filename = obs_str + "_" + filenamepart
                    fullfilename = os.path.join(datapath, filename)
                    self.assertTrue(os.path.exists(fullfilename))
                    self.assertTrue(os.path.isfile(fullfilename))

    def test_sanitizedstring(self):
        """Test sanitized string"""
        self.assertEqual(sanitizeString("L s s T"), "LSST")

    def test_validateobsname(self):
        """Test validate observatory name"""
        self.assertEqual(validateObsName("r u bin"), "LSST")
        self.assertEqual(validateObsName("rubin observatory"), "LSST")
        self.assertEqual(validateObsName("A u x Tel"), "LSST")
        self.assertEqual(validateObsName("cerro tololo"), "CTIO")
        self.assertEqual(validateObsName("pic du midi"), "PDM")
        self.assertEqual(validateObsName("Observatoire du Pic du Midi"), "PDM")
        self.assertEqual(validateObsName("mauna kea"), "OMK")
        self.assertEqual(validateObsName("mauna kea observatory"), "OMK")
        self.assertEqual(validateObsName("sea level"), "OSL")
        self.assertEqual(validateObsName("sea level observatory"), "OSL")
        self.assertIsNone(validateObsName("AaBb"))

    def test_is_obssite(self):
        """Test if observation site exist"""
        self.assertTrue(is_obssite("LSST"))
        self.assertTrue(is_obssite("a u x tel"))
        self.assertTrue(is_obssite("AUXTEL"))
        self.assertTrue(is_obssite("CTIO"))
        self.assertTrue(is_obssite("OMK"))
        self.assertFalse(is_obssite("AMissingObs"))
        self.assertTrue(is_obssite("L s s t"))

    def test_ObsAtmo(self):
        """Test ObsAtmo emulator"""
        e = ObsAtmo("LSST")
        array = e.GetAllTransparencies(
            [300.0, 400.0, 600.0, 800.0, 1000.0], 1.0, 0.0, 0.0
        )
        self.assertTrue(
            np.allclose(
                array, [0.41483073, 0.76797745, 0.95170851, 0.98476396, 0.9937721]
            ),
            msg="ObsAtmo.GetAllTransparencies error",
        )

        e = ObsAtmo("AUXTEL")
        array = e.GetAllTransparencies(
            [300.0, 400.0, 600.0, 800.0, 1000.0], 1.0, 0.0, 0.0
        )
        self.assertTrue(
            np.allclose(
                array, [0.41483073, 0.76797745, 0.95170851, 0.98476396, 0.9937721]
            ),
            msg="ObsAtmo.GetAllTransparencies error",
        )

    def test_ObsAtmoPressure(self):
        """Test ObsAtmoPressure emulator"""
        e = ObsAtmoPressure("LSST")
        array = e.GetAllTransparencies(
            [300.0, 400.0, 600.0, 800.0, 1000.0], 1.0, 0.0, 0.0
        )
        self.assertTrue(
            np.allclose(
                array, [0.41483073, 0.76797745, 0.95170851, 0.98476396, 0.9937721]
            ),
            msg="ObsAtmo.GetAllTransparencies error",
        )

        e = ObsAtmoPressure("AUXTEL")
        array = e.GetAllTransparencies(
            [300.0, 400.0, 600.0, 800.0, 1000.0], 1.0, 0.0, 0.0
        )
        self.assertTrue(
            np.allclose(
                array, [0.41483073, 0.76797745, 0.95170851, 0.98476396, 0.9937721]
            ),
            msg="ObsAtmo.GetAllTransparencies error",
        )

    def test_ObsAtmoGrid(self):
        """Test ObsAtmoGrid emulator"""
        e = ObsAtmoGrid("LSST")
        array = e.GetAllTransparencies(
            [300.0, 400.0, 600.0, 800.0, 1000.0], 1.0, 0.0, 0.0
        )
        self.assertTrue(
            np.allclose(
                array, [0.41483073, 0.76797745, 0.95170851, 0.98476396, 0.9937721]
            ),
            msg="ObsAtmo.GetAllTransparencies error",
        )

        e = ObsAtmoGrid("AUXTEL")
        array = e.GetAllTransparencies(
            [300.0, 400.0, 600.0, 800.0, 1000.0], 1.0, 0.0, 0.0
        )
        self.assertTrue(
            np.allclose(
                array, [0.41483073, 0.76797745, 0.95170851, 0.98476396, 0.9937721]
            ),
            msg="ObsAtmo.GetAllTransparencies error",
        )


# -------------------------------------------------------------------
# Special Tests for convert_dict_to_json and convert_json_to_dict
# -------------------------------------------------------------------


class TestConvertFunctions(unittest.TestCase):
    """Tests for convert_dict_to_json and convert_json_to_dict"""

    def test_convert_dict_to_json_with_numpy_array(self):
        data = {"array": np.array([1, 2, 3]), "value": 42}
        result = convert_dict_to_json(data)
        self.assertEqual(result, {"array": [1, 2, 3], "value": 42})

    def test_convert_dict_to_json_with_scalars(self):
        data = {"int": 5, "float": 3.14, "str": "hello", "bool": True, "none": None}
        result = convert_dict_to_json(data)
        self.assertEqual(result, data)

    def test_convert_dict_to_json_invalid_type(self):
        data = {"unsupported": set([1, 2, 3])}
        with self.assertRaises(TypeError):
            convert_dict_to_json(data)

    def test_convert_dict_to_json_invalid_input(self):
        with self.assertRaises(TypeError):
            convert_dict_to_json(["not", "a", "dict"])

    def test_convert_json_to_dict_with_list(self):
        data_json = {"array": [1, 2, 3], "value": 42}
        result = convert_json_to_dict(data_json)
        self.assertTrue(np.array_equal(result["array"], np.array([1, 2, 3])))
        self.assertEqual(result["value"], 42)

    def test_convert_json_to_dict_with_nonlist_values(self):
        data_json = {"a": 1, "b": "test", "c": None}
        result = convert_json_to_dict(data_json)
        self.assertEqual(result, data_json)

    def test_convert_json_to_dict_invalid_input(self):
        with self.assertRaises(TypeError):
            convert_json_to_dict("not a dict")

    def test_round_trip_conversion(self):
        """Test that converting dict -> json -> dict preserves data"""
        original = {"arr": np.arange(4), "val": 7.5, "txt": "ok"}
        json_data = convert_dict_to_json(original)
        recovered = convert_json_to_dict(json_data)
        self.assertTrue(np.array_equal(recovered["arr"], original["arr"]))
        self.assertEqual(recovered["val"], original["val"])
        self.assertEqual(recovered["txt"], original["txt"])


if __name__ == "__main__":
    unittest.main()
