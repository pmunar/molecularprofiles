"""
Calculates the properties of moist air, including:
* density
* compressibility
* enhancement factor
* saturation vapor pressure
* molar fraction from relative humidity
* partial pressure from relative humidity

Adapted from MHumidity, written by Markus Gaug <markus.gaug@uab.cat>, 04/2013
"""

import numpy as np
from math import exp, sqrt, log10

MOLAR_MASS_WATER_VAPOR = 0.018015   # molar mass of water vapor [kg/mol]
GAS_CONSTANT = 8.31451              # gas constant [J/mol/K]


def Compressibility(p, T, Xw):
    """
    Calculates the compressibility of moist air, according to:
    R.S. Davis, "Equation for the determination of the density of moist air"
    Metrologia, 29 (1992) 67-70

    See also Eq. 16 in:
    C. Tomasi, V. Vitale, B. Petkov, A. Lupi, A. Cacciari
    "Improved algorithm for calculations of Rayleigh-scattering optical depth
    in standard atmospheres", Applied Optics 44 Nr. 16 (2005) 3320

    Args:
        p  : pressure in hPa
        T  : temperature in K
        Xw : molar fraction of water vapor
    Returns:
        compressibility of moist air (dimensionless constant, 0 < Z < 1)
    """
    Tc = T - 273.15   # temperature in deg. C
    pT = 100*p/T      # ratio of pressure to temperature, with pressure in Pascals
    a0 = 1.58123E-6   # K Pa^-1
    a1 = -2.9331E-8   # Pa^-1
    a2 = 1.1043E-10   # K^-1 Pa^-1
    b0 = 5.707E-6     # K Pa^-1
    b1 = -2.051E-8    # Pa^-1
    c0 = 1.9898E-4    # K Pa^-1
    c1 = -2.376E-6    # Pa^-1
    d0 = 1.83E-11     # K^2 Pa^-2
    d1 = -7.65E-9     # K^2 Pa^-2
    Z = 1 - pT*(a0 + a1*Tc + a2*Tc*Tc + b0*Xw + b1*Xw*Tc + c0*Xw*Xw + c1*Xw*Xw*Tc) + pT*pT*(d0 + d1*Xw*Xw)
    return Z


def EnhancementFactor(p, T):
    """
    Calculates the enhancement factor of water vapor in air.

    Calculated according to Eq. 14 of:
    C. Tomasi, V. Vitale, B. Petkov, A. Lupi, A. Cacciari
    "Improved algorithm for calculations of Rayleigh-scattering optical depth
    in standard atmospheres", Applied Optics 44 Nr. 16 (2005) 3320

    Args:
        p : pressure in Pa
        T : temperature in K
    Returns:
        enhancement factor (dimensionless constant)
    """
    Tc = T - 273.15   # temparture in deg. C
    f = 1.00062 + 3.14E-8*p + 5.6E-7*pow(Tc, 2)
    return f


# TODO: This function needs a simple unit test to check which function is called
def SaturationVaporPressure(T):
    if T > 273.15:
        return SaturationVaporPressureDavis(T)
    else:
        return SaturationVaporPressureGoffGratch(T)


def SaturationVaporPressureDavis(T):
    """
    Calculates the vapor pressure at saturation, according to:
    R.S. Davis, "Equation for the determination of the density of moist air"
    Metrologia, 29 (1992) 67-70

    See also Eq. 15 in:
    C. Tomasi, V. Vitale, B. Petkov, A. Lupi, A. Cacciari
    "Improved algorithm for calculations of Rayleigh-scattering optical depth
    in standard atmospheres", Applied Optics 44 Nr. 16 (2005) 3320

    Args:
        T : temperature in K
    Returns:
        saturation vapor pressure in hPa
    """
    E = np.exp(1.2378847e-5*T*T - 1.9121316e-2*T + 33.93711047 - 6343.1645/T)
    return E/100   # divide by 100 for hPa


# TODO: Goff-Gratch doesn't seem to be as good as Buck, or maybe Wexler
def SaturationVaporPressureGoffGratch(T):
    """
    Calculates the vapor pressure at saturation, according to:
    Smithsonian Tables (1984); after Goff and Gratch (1946).
    See here: http://cires.colorado.edu/~voemel/vp.html

    This equation is recommended for temperatures below 0 deg. C.

    Args:
        T : temperature in K
    Returns:
        saturation vapor pressure in hPa
    """
    theta = 373.16/T   # ratio of steam point (100 deg C) to temperature
    c = [-7.90298*(theta - 1),
         5.02808*log10(theta),
         -1.3816e-7*(pow(10, 11.344*(1 - 1/theta)) - 1),
         8.1328e-3*(pow(10, -3.49149*(theta - 1)) - 1),
         log10(1013.246)]
    log10_ew = np.sum(c)
    ew = pow(10, log10_ew)
    return ew


# TODO: Check this function and write unit test
def SaturationVaporPressureOverWater(T):
    """
    Calculates the vapor pressure at saturation over water, according to IAPWS:
    International Association for the Properties of Water and Steam,
    Peter H. Huang, "New equations for water vapor pressure in the temperature
    range -100 deg. C to 100 deg. C for use with the 1997 NIST/ASME steam tables"
    Papers and abstracts from the third international symposium on humidity and
    moisture, Vol. 1, p. 69-76, National Physical Laboratory, Teddington,
    Middlesex, UK, April 1998.
    See also: http://cires.colorado.edu/~voemel/vp.html

    Args:
        T : temperature in K
    Returns:
        saturation vapor pressure in hPa
    """
    omega = T - 2.38555575678E-01/(T - 6.50175348448E+02)
    omega2 = omega*omega
    A =                    omega2 + 1.16705214528E+03*omega - 7.24213167032E+05
    B = -1.70738469401E+01*omega2 + 1.20208247025E+04*omega - 3.23255503223E+06
    C =  1.49151086135E+01*omega2 - 4.82326573616E+03*omega + 4.05113405421E+05
    X = -B + sqrt(B*B - 4*A*C)
    return 1e4*pow(2*C/X, 4)


# TODO: Check this function and write unit test
def SaturationVaporPressureOverIce(T):
    """
    Calculates the vapor pressure at saturation over ice, according to IAPWS:
    International Association for the Properties of Water and Steam,
    Peter H. Huang, "New equations for water vapor pressure in the temperature
    range -100 deg. C to 100 deg. C for use with the 1997 NIST/ASME steam tables"
    Papers and abstracts from the third international symposium on humidity and
    moisture, Vol. 1, p. 69-76, National Physical Laboratory, Teddington,
    Middlesex, UK, April 1998.
    See also: http://cires.colorado.edu/~voemel/vp.html

    Args:
        T : temperature in K
    Returns:
        saturation vapor pressure in hPa
    """
    theta = T/273.16
    Y = -13.928169*(1 - pow(theta, -1.5)) + 34.7078238*(1 - pow(theta, -1.25))
    return 6.11657*exp(Y)


def MolarFractionWaterVapor(p, T, RH):
    """
    Calculates the molar fraction of water vapor in moist air.

    See the text above Eq. 14 of:
    C. Tomasi, V. Vitale, B. Petkov, A. Lupi, A. Cacciari
    "Improved algorithm for calculations of Rayleigh-scattering optical depth
    in standard atmospheres", Applied Optics 44 Nr. 16 (2005) 3320

    Args:
        p  : pressure in hPa
        T  : temperature in K
        RH : relative humidity in percent
    Returns:
        Xw : molar fraction of water vapor in moist air (dimensionless)
    """
    f = EnhancementFactor(100*p, T)
    psv = SaturationVaporPressure(T)
    Xw = f * RH/100 * psv/p   # Xw = f * h * E(T)/p
    return Xw


def DensityMoistAir(p, T, Z, Xw, C):
    """
    Density equation of moist air, according to:
    R.S. Davis, "Equation for the determination of the density of moist air"
    Metrologia, 29 (1992) 67-70

    Args:
        p  : pressure in hPa (beware, different unit than in Davis!)
        T  : temperature in Kelvin
        Z  : compressibility (see Compressibility() in this module)
        Xw : molar fraction of water vapor
        C  : CO2 volume concentration in ppmv (different unit than in Davis!)
    Returns:
        density of moist air (kg m^-3)
    """
    p *= 100   # convert pressure to Pa
    R = GAS_CONSTANT
    Mw = MOLAR_MASS_WATER_VAPOR
    Ma = 1e-3*(28.9635 + 12.011e-6*(C - 400))   # molar mass of dry air [kg/mol]
    rho = p*Ma/(Z*R*T) * (1 - Xw*(1 - Mw/Ma))   # Tomasi eq. 12
    return rho


# TODO: This function has no unit test
def PartialPressureWaterVapor(T, RH):
    """
    Calculates the partial pressure of water vapor in the air.

    Args:
        T  : temperature in K
        RH : relative humidity in percent
    Returns:
        ew : water vapor partial pressure in hPa
    """
    # water vapor pressure: e = h * E(T)
    ew = (RH/100)*SaturationVaporPressureDavis(T)
    return ew
