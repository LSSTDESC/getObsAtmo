"""Initialisation du package getObsAtmo."""

from ._version import __version__

from .getObsAtmo import (
    Dict_Of_sitesAltitudes,
    Dict_Of_sitesPressures,
    Dict_Of_sitesMoleculesNDensity,
    Dict_Of_sitesAirWeight,
    Dict_Of_sitesAliases,
    file_data_dict,
    _getPackageDir,
    get_obssite_keys,
    sanitizeString,
    validateObsName,
    is_obssite,
    ObsAtmo,
    ObsAtmoPressure,
    ObsAtmoGrid,
    convert_dict_to_json,
    convert_json_to_dict,
)

__all__ = [
    "__version__",
    "Dict_Of_sitesAltitudes",
    "Dict_Of_sitesPressures",
    "Dict_Of_sitesMoleculesNDensity",
    "Dict_Of_sitesAirWeight",
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
    "convert_dict_to_json",
    "convert_json_to_dict",
]
