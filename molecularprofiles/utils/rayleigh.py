"""
A small Rayleigh-scattering program based on:
    C. Tomasi, V. Vitale, B. Petkov, A. Lupi, A. Cacciari
    "Improved algorithm for calculations of Rayleigh-scattering optical depth in standard atmospheres"
    Applied Optics 44 Nr. 16 (2005) 3320

The calculation of refractive index is based on:
    P.E. Ciddor, "Refractive index of air: new equations for hte visible and near infrared",
    Applied Optics 35 (1996) 1566

    P.E. Ciddor, "Refractive index of air: 3. The roles of CO2, H20 and refractivity virals",
    Applied Optics 41 (2002) 2292

The principal King factor formula is based on:
    D.R. Bates, "Rayleigh scattering by air",
    Planet. Space Sci. 32 (1984) 785

    B.A. Bodhaine, N.B. Wood, E.G. Dutton, J.R. Slusser, "On Rayleigh optical depth calculations",
    J. Atmosph. Osceanic Technol. 16 (1999) 1854

The calculation of the Chandrasekhar phase function is based on:
    S. Chandrasekhar, Radiative Transfer, Dover Publications, 1960.

    E.J. McCartney, "Optics of the Atmosphere. Scattering by Molecules and Particles"
    Wiley & Sons, New York, 1977.

Adapted from MRayleigh, written by Markus Gaug <markus.gaug@uab.cat>, 04/2013
.. moduleauthor:: Scott Griffiths <sgriffiths@ifae.es>
"""

from __future__ import division, print_function
import matplotlib.pyplot as plt
from math import pi, cos, log
import numpy as np
import molecularprofiles.utils.humidity as humidity


class Rayleigh:
    # Class Constants
    Ns = 2.546899e19   # [cm^-3] molecular number density for standard air conditions
    ps = 1013.25       # [hPa]   standard air pressure
    Ts = 288.15        # [K]     standard air temperature
    RHs = 45.9         # [%]     standard air rel. humidity
    Cs = 385           # [mmpv]  standard air CO2 concentration

    # TODO: do we need these?
    # p_mean = 788.2                      # [hPa] mean air pressure at MAGIC site
    # T_mean = 273.15 + 8.89              # [K]   mean air temperature at MAGIC site
    # LIDAR_ratio_mol = 8.37758           # for molecules: 8pi/3

    def __init__(self, wavelength, p=ps, T=Ts, RH=RHs, C=Cs):
        """
        Constructor for Rayleigh-scattering class.

        Args:
            wavelength : wavelength of light [nm]
            p  : pressure [hPa]
            T  : temperature [K]
            RH : relative humidity [%]
            C  : CO2 concentration [ppmv]
        """
        # check inputs for bad values
        if wavelength < 200 or wavelength > 4000:
            raise ValueError("Wavelength range only from 200 nm - 4 micrometer allowed.")
        # if p < 100 or p > 1400:
        #     raise ValueError("Pressure only in range 100 - 1400 hPa allowed.")
        if p < 0 or p > 1400:   # changed to make compatible with gdas_entry
            raise ValueError("Pressure only in range 0 - 1400 hPa allowed.")
        if T < 200 or T > 373.15:
            raise ValueError("Temperatures only in range 200 - 373 K allowed.")
        if RH < 0 or RH > 100:
            raise ValueError("Relative humity must lie between 0 - 100.")
        if C < 200 or C > 1000:
            raise ValueError("CO2 concentrations only in range 200 - 1000 ppmv allowed.")

        self.wavelength = wavelength   # [nm]    wavelenght of light
        self.p = p                     # [hPa]   air pressure
        self.T = T                     # [K]     air temperature
        self.RH = RH                   # [%]     relative humidity of air
        self.C = C                     # [ppmv]  CO2 concentration of air

        self.calculate_sigma()   # [cm^-2] single molecule scattering coefficient
        self.calculate_beta()    # [km^-1] volume scattering coefficient

    def calculate_sigma(self):
        """
        Calculate the total Rayleigh scattering cross section per moleculem, sigma,
        in units of cm^-2.
        """
        wavelength_cm = self.wavelength*1e-7
        l4 = pow(wavelength_cm, 4)

        # calculate the refractive index of air
        n = self.calculate_n()
        n2 = n*n

        # calculate molecular number density in units of cm^-3 (Tomasi eq. 3)
        N = self.Ns * self.p/self.ps * self.Ts/self.T
        N2 = N*N

        pi3 = pow(pi, 3)
        F = self.king()
        self.sigma = 24*pi3*(n2 - 1)**2 / (l4*N2*(n2 + 2)**2) * F   # units of cm^-2
        return self.sigma

    def calculate_beta(self):
        """
        Calculates the monochromatic volume coefficient for the total molecular
        scattering in cloudless air, beta, in units of km^-1
        """
        # calculate molecular number density in units of cm^-3 (Tomasi eq. 3)
        N = self.Ns * self.p/self.ps * self.Ts/self.T

        try:   # check that sigma has been calculated
            self.sigma
        except NameError:
            self.calculate_sigma()

        # calculate beta in units of km^-1 (Tomasi eq. 2)
        self.beta = 1e5*N*self.sigma
        return self.beta

    def calculate_n(self):
        """
        Ciddor formula for calculation of refractive index in moist air.
        The obtained refractive index is precise to 1e-7.

        Cross-checked with:
        http://emtoolbox.nist.gov/Wavelength/Documentation.asp#IndexofRefractionofAir

        Returns:
            n : index of refraction of moist air
        """
        l = pow(self.wavelength/1000, -2)   # convert wavelength to micrometers

        # refractive index of standard dry air (e = 0) according to Ciddor, with 450 ppmv CO2
        n_as = 1e-8*(5792105/(238.0185 - l) + 167917/(57.362 - l)) + 1   # Tomasi eq. 17

        # refractive index of dry air at standard p and T, for given C (e = 0)
        n_axs = (1 + 0.534e-6*(self.C - 450))*(n_as - 1) + 1   # Tomasi eq. 18

        # refractive index of pure water vapor at standard T and e (T* = 293.15 K = 20 C, and e* = 1333 Pa)
        n_ws = 1.022e-8*(295.235 + 2.6422*l - 0.032380*l*l + 0.004028*l*l*l) + 1   # Tomasi eq. 19

        # calculate the respective densities (see Tomasi et al., pp. 3325 ff)
        R = 8.31451                                          # gas constant [J/mol/K]
        Ma = 1e-3*(28.9635 + 12.011e-6*(self.C - 400))       # molar mass of dry air [kg/mol]
        Mw = 0.018015                                        # molar mass of water vapor [kg/mol]
        Xw = humidity.MolarFractionWaterVapor(self.p, self.T, self.RH)   # molar fraction of water vapor in moist air
        Za = humidity.Compressibility(self.ps, self.Ts, 0)   # compressibility of dry air
        Zw = humidity.Compressibility(13.33, 293.15, 1)      # compressibility of pure water vapor
        Zm = humidity.Compressibility(self.p, self.T, Xw)    # compressibility of moist air

        # density of dry air at standard p and T
        rho_axs = humidity.DensityMoistAir(self.ps, self.Ts, Za, 0, self.C)

        # density of pure water vapor at at standard T and e (T* = 293.15 K = 20 C, and e* = 1333 Pa)
        rho_ws = humidity.DensityMoistAir(13.33, 293.15, Zw, 1, self.C)

        # density of the dry component of moist air
        rho_a = (100*self.p)*Ma*(1 - Xw)/(Zm*R*self.T)

        # density of the water vapor component of moist air
        rho_w = (100*self.p)*Mw*Xw/(Zm*R*self.T)

        n = 1 + (rho_a/rho_axs)*(n_axs - 1) + (rho_w/rho_ws)*(n_ws - 1)   # Ciddor eq. 5, Tomasi eq. 11
        return n

    def king(self):
        """
        Calculates the current best estimate of the King factor of moist air.

        The King factor is used to take into account effects due to the anisotropic
        properties of air molecules since anisotropic molecules scatter more radiation
        at 90 degrees scattering angles than isotropic molecules with the same index
        of refraction.

        Precision not stated in Tomasi et al., but probably better than 1e-4.
        Effects of RH are of the order of several times 1e-4.

        Returns:
            F : King factor [dimensionless]
        """
        l = pow(self.wavelength/1000, -2)   # convert to micrometers
        e = humidity.PartialPressureWaterVapor(self.T, self.RH)   # water vapor partial pressure [hPa]

        F1 = 1.034 + 3.17e-4*l                   # partial King factor for N2 molecules
        F2 = 1.096 + 1.385e-3*l + 1.448e-4*l*l   # partial King factor for O2 molecules
        F3 = 1.00                                # partial King factor for Ar molecules
        F4 = 1.15                                # partial King factor for CO2 molecules
        F5 = 1.001                               # partial King factor for water vapor

        c1 = 0.78084       # N2
        c2 = 0.20946       # O2
        c3 = 0.00934       # Ar
        c4 = 1e-6*self.C   # CO2
        c5 = e/self.p      # water vapor mixing ratio

        F = (c1*F1 + c2*F2 + c3*F3 + c4*F4 + c5*F5)/(c1 + c2 + c3 + c4 + c5)   # Tomasi eq. 22
        return F

    def phase_function(self, angle):
        """
        Calculates the Chandrasekhar phase function.

        Args:
            angle : scattering angle in radians

        Returns:
            Chandrasekhar phase function for scattering of natural light
        """
        rho = self.depolarization()

        # need to solve Chandrasekhar eq. 254 for gamma as a function of rho
        f1 = (2 + 2*rho)/(2 + rho)
        f2 = (1 - rho)/(1 + rho)
        return 0.75*f1*(1 + f2*pow(cos(angle), 2))   # Chandrasekhar eq. 255

    def depolarization(self):
        """
        Current best estimate of the depolarization factor of moist air.

        Precision not stated in Tomasi et al., but probably better than 1e-4.
        Effects of RH are of the order of several times 1e-4.
        """
        F = self.king()
        return 6*(F - 1)/(3 + 7*F)   # Tomasi eq. 5, solved for rho

    # TODO: where does this come from? Needs unit test
    def dbeta_domega(self, angle):
        """
        Get the back-scattering coefficient for a given scattering angle.

        Args:
            angle : scattering angle in radians

        Returns:
            back-scattering coefficient in units of km^-1
        """
        try:   # check that beta has been calculated
            self.beta
        except NameError:
            self.calculate_beta()

        phase_func = self.phase_function(angle)
        return phase_func*self.beta/(4*pi)

    def print_params(self):
        """Prints Rayleigh scattering parameters."""
        print("Wavelength:              %.1f\t nm" % self.wavelength)
        print("Air Pressure:            %.2f\t hPa" % self.p)
        print("Air Temperature:         %.2f\t K" % self.T)
        print("Rel. Humidity:           %.1f\t %%" % self.RH)
        print("CO2 concentration:       %.0f\t ppmv" % self.C)
        print("Refractive Index:        %.8f" % self.calculate_n())
        print("King Factor:             %.5f" % self.king())
        print("Depolarization:          %.4f" % self.depolarization())
        print("Mol. cross section:      %.6g\t cm^-2" % self.sigma)
        print("Volume scattering coeff: %.5g\t km^-1" % self.beta)

    def plot(self, param, nbins=20):
        """
        Convenience function to plot one of the following parameters:
        wavelength, phase, p, T, RH, C
        """
        plot_dict = {'wavelength': self.plot_wavelength,
                     'phase': self.plot_phase_function,
                     'p': self.plot_pressure,
                     'T': self.plot_temperature,
                     'RH': self.plot_rel_humidity,
                     'C': self.plot_co2_concentration}
        plotter = plot_dict[param]
        plotter(nbins)

    def plot_wavelength(self, nbins=20):
        """Plots sigma and beta vs. wavelength."""
        wavelengths = 200 + np.arange(nbins)/nbins*(4000 - 200)
        y_sigma = []; y_beta = []
        for l in wavelengths:
            r = Rayleigh(l, self.p, self.T, self.RH, self.C)
            y_sigma.append(r.sigma)
            y_beta.append(r.beta)

        p_str = r'Pressure p={:.2f} hPa'.format(self.p)
        T_str = r'Temperature T={:.2f} K'.format(self.T)
        RH_str = r'Rel. Humidity RH={:.1f}\%'.format(self.RH)
        CO2_str = r'CO$_2$ Concentration={:.0f} ppmv'.format(self.C)
        box_text = '\n'.join([p_str, T_str, RH_str, CO2_str])

        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

        axes[0].semilogy(wavelengths, y_sigma, 'ro')
        axes[0].set_xlabel(r"$\lambda$ (nm)")
        axes[0].set_ylabel(r"$\sigma$ (cm$^{-2}$)")
        axes[0].set_title("Total Scattering Coefficient vs. Wavelength")
        axes[0].text(0.45, 0.89, box_text, transform=axes[0].transAxes, verticalalignment='top', linespacing=2)
        axes[0].minorticks_on()

        axes[1].semilogy(wavelengths, y_beta, 'bs')
        axes[1].set_xlabel(r"$\lambda$ (nm)")
        axes[1].set_ylabel(r"$\beta$ (km$^{-1}$)")
        axes[1].set_title("Volume Scattering Coefficient vs. Wavelength")
        axes[1].text(0.45, 0.89, box_text, transform=axes[1].transAxes, verticalalignment='top', linespacing=2)
        axes[1].minorticks_on()

        fig.subplots_adjust(left=0.07, right=0.96, bottom=0.09, top=0.93)

    def plot_phase_function(self, nbins=20):
        """Plots the phase function and beta(theta) vs. angle."""
        angles = np.arange(nbins)/(nbins-1)*180
        y_sigma = []; y_beta = []
        for angle in angles:
            f = self.phase_function(np.deg2rad(angle))
            y_sigma.append(f)
            y_beta.append(f*self.beta/(4*pi))

        l_str = r'Wavelength $\lambda$={:.1f} nm'.format(self.wavelength)
        p_str = r'Pressure p={:.2f} hPa'.format(self.p)
        T_str = r'Temperature T={:.2f} K'.format(self.T)
        RH_str = r'Rel. Humidity RH={:.1f}\%'.format(self.RH)
        CO2_str = r'CO$_2$ Concentration={:.0f} ppmv'.format(self.C)
        box_text = '\n'.join([l_str, p_str, T_str, RH_str, CO2_str])

        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

        axes[0].plot(angles, y_sigma, 'ro')
        axes[0].set_xlabel(r"$\theta$ (deg)")
        axes[0].set_ylabel(r"$P(\theta)$")
        axes[0].set_title("Phase Function")
        axes[0].text(0.5, 0.95, box_text, transform=axes[0].transAxes, verticalalignment='top', horizontalalignment='center', linespacing=2)
        axes[0].minorticks_on()

        axes[1].plot(angles, y_beta, 'bs')
        axes[1].set_xlabel(r"$\theta$ (deg)")
        axes[1].set_ylabel(r"$\beta(\theta)$ (km$^{-1}$)")
        axes[1].set_title("Volume Scattering Coefficient vs. Scattering Angle")
        axes[1].text(0.5, 0.95, box_text, transform=axes[1].transAxes, verticalalignment='top', horizontalalignment='center', linespacing=2)
        axes[1].minorticks_on()

        fig.subplots_adjust(left=0.07, right=0.96, bottom=0.09, top=0.93)

    def plot_pressure(self, nbins=20):
        """Plots sigma and beta vs. pressure."""
        pressures = 100 + np.arange(nbins)/nbins*(1025 - 100)
        y_sigma = []; y_beta = []
        for p in pressures:
            r = Rayleigh(self.wavelength, p, self.T, self.RH, self.C)
            y_sigma.append(r.sigma)
            y_beta.append(r.beta)

        l_str = r'Wavelength $\lambda$={:.1f} nm'.format(self.wavelength)
        T_str = r'Temperature T={:.2f} K'.format(self.T)
        RH_str = r'Rel. Humidity RH={:.1f}\%'.format(self.RH)
        CO2_str = r'CO$_2$ Concentration={:.0f} ppmv'.format(self.C)
        box_text = '\n'.join([l_str, T_str, RH_str, CO2_str])

        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

        axes[0].plot(pressures, y_sigma, 'ro')
        axes[0].set_xlabel(r"P (hPa)")
        axes[0].set_ylabel(r"$\sigma$ (cm$^{-2}$)")
        axes[0].set_title("Total Scattering Coefficient vs. Pressure")
        axes[0].text(0.7, 0.25, box_text, transform=axes[0].transAxes, verticalalignment='top', horizontalalignment='center', linespacing=2)
        axes[0].minorticks_on()

        axes[1].plot(pressures, y_beta, 'bs')
        axes[1].set_xlabel(r"P (hPa)")
        axes[1].set_ylabel(r"$\beta$ (km$^{-1}$)")
        axes[1].set_title("Volume Scattering Coefficient vs. Pressure")
        axes[1].text(0.7, 0.25, box_text, transform=axes[1].transAxes, verticalalignment='top', horizontalalignment='center', linespacing=2)
        axes[1].minorticks_on()

        fig.subplots_adjust(left=0.07, right=0.96, bottom=0.09, top=0.93)

    def plot_temperature(self, nbins=20):
        """Plots sigma and beta vs. temperature."""
        temperatures = 200 + np.arange(nbins)/nbins*(380 - 200)
        y_sigma = []; y_beta = []
        for T in temperatures:
            r = Rayleigh(self.wavelength, self.p, T, self.RH, self.C)
            y_sigma.append(r.sigma)
            y_beta.append(r.beta)

        l_str = r'Wavelength $\lambda$={:.1f} nm'.format(self.wavelength)
        p_str = r'Pressure p={:.2f} hPa'.format(self.p)
        RH_str = r'Rel. Humidity RH={:.1f}\%'.format(self.RH)
        CO2_str = r'CO$_2$ Concentration={:.0f} ppmv'.format(self.C)
        box_text = '\n'.join([l_str, p_str, RH_str, CO2_str])

        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

        axes[0].plot(temperatures, y_sigma, 'ro')
        axes[0].set_xlabel(r"T (K)")
        axes[0].set_ylabel(r"$\sigma$ (cm$^{-2}$)")
        axes[0].set_title("Total Scattering Coefficient vs. Temp.")
        axes[0].text(0.3, 0.25, box_text, transform=axes[0].transAxes, verticalalignment='top', horizontalalignment='center', linespacing=2)
        axes[0].minorticks_on()

        axes[1].plot(temperatures, y_beta, 'bs')
        axes[1].set_xlabel(r"T (K)")
        axes[1].set_ylabel(r"$\beta$ (km$^{-1}$)")
        axes[1].set_title("Volume Scattering Coefficient vs. Temp.")
        axes[1].text(0.3, 0.25, box_text, transform=axes[1].transAxes, verticalalignment='top', horizontalalignment='center', linespacing=2)
        axes[1].minorticks_on()

        fig.subplots_adjust(left=0.07, right=0.96, bottom=0.09, top=0.93)

    def plot_rel_humidity(self, nbins=20):
        """Plots sigma and beta vs. relative humidity."""
        rel_humidities = np.arange(nbins)/(nbins-1)*100
        y_sigma = []; y_beta = []
        for RH in rel_humidities:
            r = Rayleigh(self.wavelength, self.p, self.T, RH, self.C)
            y_sigma.append(r.sigma)
            y_beta.append(r.beta)

        l_str = r'Wavelength $\lambda$={:.1f} nm'.format(self.wavelength)
        p_str = r'Pressure p={:.2f} hPa'.format(self.p)
        T_str = r'Temperature T={:.2f} K'.format(self.T)
        CO2_str = r'CO$_2$ Concentration={:.0f} ppmv'.format(self.C)
        box_text = '\n'.join([l_str, p_str, T_str, CO2_str])

        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

        axes[0].plot(rel_humidities, y_sigma, 'ro')
        axes[0].set_xlabel(r"RH (\%)")
        axes[0].set_ylabel(r"$\sigma$ (cm$^{-2}$)")
        axes[0].set_title("Total Scattering Coeff. vs. Rel. Humidity")
        axes[0].text(0.7, 0.95, box_text, transform=axes[0].transAxes, verticalalignment='top', horizontalalignment='center', linespacing=2)
        axes[0].minorticks_on()

        axes[1].plot(rel_humidities, y_beta, 'bs')
        axes[1].set_xlabel(r"RH (\%)")
        axes[1].set_ylabel(r"$\beta$ (km$^{-1}$)")
        axes[1].set_title("Volume Scattering Coeff. vs. Rel. Humidity")
        axes[1].text(0.7, 0.95, box_text, transform=axes[1].transAxes, verticalalignment='top', horizontalalignment='center', linespacing=2)
        axes[1].get_yaxis().get_major_formatter().set_useOffset(False)
        axes[1].minorticks_on()

        fig.subplots_adjust(left=0.07, right=0.96, bottom=0.09, top=0.93)

    def plot_co2_concentration(self, nbins=20):
        """Plots sigma and beta vs. C02 concentration."""
        concentrations = 300 + np.arange(nbins)/(nbins-1)*(450 - 300)
        y_sigma = []; y_beta = []
        for C in concentrations:
            r = Rayleigh(self.wavelength, self.p, self.T, self.RH, C)
            y_sigma.append(r.sigma)
            y_beta.append(r.beta)

        l_str = r'Wavelength $\lambda$={:.1f} nm'.format(self.wavelength)
        p_str = r'Pressure p={:.2f} hPa'.format(self.p)
        T_str = r'Temperature T={:.2f} K'.format(self.T)
        RH_str = r'Rel. Humidity RH={:.1f}\%'.format(self.RH)
        box_text = '\n'.join([l_str, p_str, T_str, RH_str])

        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

        axes[0].plot(concentrations, y_sigma, 'ro')
        axes[0].set_xlabel(r"CO$_2$ (ppmv)")
        axes[0].set_ylabel(r"$\sigma$ (cm$^{-2}$)")
        axes[0].set_title(r"Total Scattering Coeff. vs. CO$_2$ Concent.")
        axes[0].text(0.7, 0.25, box_text, transform=axes[0].transAxes, verticalalignment='top', horizontalalignment='center', linespacing=2)
        axes[0].ticklabel_format(axis='y', style='sci', useOffset=False, scilimits=(0, 0))
        axes[0].minorticks_on()

        axes[1].plot(concentrations, y_beta, 'bs')
        axes[1].set_xlabel(r"CO$_2$ (ppmv)")
        axes[1].set_ylabel(r"$\beta$ (km$^{-1}$)")
        axes[1].set_title(r"Volume Scattering Coeff. vs. CO$_2$ Concent.")
        axes[1].text(0.7, 0.25, box_text, transform=axes[1].transAxes, verticalalignment='top', horizontalalignment='center', linespacing=2)
        axes[1].ticklabel_format(axis='y', style='sci', useOffset=False, scilimits=(0, 0))
        axes[1].minorticks_on()

        fig.subplots_adjust(left=0.07, right=0.96, bottom=0.09, top=0.93)
