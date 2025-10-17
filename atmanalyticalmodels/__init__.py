"""
Module `atmanalyticalmodels`
====================

This package contains analytical models used to validate or cross-check
the `getObsAtom` atmospheric transmission computations.

It includes:
- Temperature and pressure as a function of altitude
- Mass and particle density profiles
- Rayleigh and aerosol optical depth formulations

Typical import:
    >>> from atmanalyticalmodels import (
    ...     Temperature_adiabatic,
    ...     Pressure_isothermal,
    ...     RayOptDepth_adiabatic
    ... )

Or:
    >>> import atmanalyticalmodels as atm
    >>> atm.Pressure_isothermal(2000)
"""

from .libatmscattering import (
    # Physical constants
    N_A,
    R,
    g0,
    M_air,
    M_air_dry,
    M_h2o,
    P0,
    T0,
    L,
    LSST_Altitude,
    CTIO_Altitude,
    OHP_Altitude,
    PDM_Altitude,
    altitude0,
    # Main functions
    Temperature_adiabatic,
    Pressure_isothermal,
    Pressure_adiabatic,
    MassDensity_isothermal,
    MassDensity_adiabatic,
    MassDensity_adiabatic_humid,
    AtomeDensity_isothermal,
    AtomeDensity_adiabatic,
    XDepth_isothermal,
    XDepth_adiabatic,
    RayOptDepth_adiabatic,
    RayOptDepth_isothermal,
    RayOptDepth2_adiabatic,
    RayOptDepth2_isothermal,
    RayOptDepthXD,
    AeroOptDepth,
)

__all__ = [
    # Physical constants
    "N_A",
    "R",
    "g0",
    "M_air",
    "M_air_dry",
    "M_h2o",
    "P0",
    "T0",
    "L",
    "LSST_Altitude",
    "CTIO_Altitude",
    "OHP_Altitude",
    "PDM_Altitude",
    "altitude0",
    # Main functions
    "Temperature_adiabatic",
    "Pressure_isothermal",
    "Pressure_adiabatic",
    "MassDensity_isothermal",
    "MassDensity_adiabatic",
    "MassDensity_adiabatic_humid",
    "AtomeDensity_isothermal",
    "AtomeDensity_adiabatic",
    "XDepth_isothermal",
    "XDepth_adiabatic",
    "RayOptDepth_adiabatic",
    "RayOptDepth_isothermal",
    "RayOptDepth2_adiabatic",
    "RayOptDepth2_isothermal",
    "RayOptDepthXD",
    "AeroOptDepth",
]
