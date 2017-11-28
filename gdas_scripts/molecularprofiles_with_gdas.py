import os.path
import sys

sys.path.append('/Volumes/Segon_HD/molecularprofiles/ecmwf/scripts_atca')
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
from scipy.interpolate import interp1d
from grib_utils import *
from readgdas import *
from readecmwf import *
from plot_settings import settings
from magic_winter_params import heightmw, rhomw
from meteorological_params import *

settings()


class EcmwfMolecularProfile:
    def __init__(self, file_ecmwf, year, file_gdas=None, tag_name='myplots', epoch_text='all',
                 label_ecmwf='ECMWF', label_gdas='GDAS', observatory='north'):

        """
        This class provides with a series of functions to analyze the quality of the ECMWF and GDAS data for both
        CTA-North and CTA-South.

        :param file_ecmwf: txt file containing the ECMWF data (string)
        :param file_gdas: txt file containing the GDAS data (string)
        :param year_text: year or years to be analyzed (just to put into some titles) (string)
        :param tag_name: name to be given to the output files (string)
        :param epoch_text: season to be analyzed. Valid values are: "all", "winter", "summer", "intermediate" (string)
        :param label_ecmwf: label to be put in some of the output plots (string)
        :param label_gdas: label to be put in some of the output plots (string)
        :param grib_file: grib file containing the original ECMWF model data (in case file_ecmwf does not exist) (string)
        :param observatory: valid options are: "north", "south"

        Functions within this class:

        retrievedata
        retrievedatafromfile: same as retrievedata but in this case it also performs the calculation of the averages
                              for several variables
        plot_average_at_15km
        plot_differences_averages
        plot_differences_wrt_magic
        plot_differences_wrt_prod3
        plot_models_comparison
        print_to_text_file

        """

        self.file_ecmwf = file_ecmwf
        self.file_gdas = file_gdas
        self.year = year
        self.epoch_text = epoch_text
        self.tag_name = tag_name
        self.label_ecmwf = label_ecmwf
        self.label_gdas = label_gdas
        self.observatory = observatory
        self.output_plot_name = self.tag_name + '_' + self.epoch_text

        # == INPUT FILE NAMES & Parameters: ==
        # file_gdas = '/Volumes/Segon_HD/gdas/GDAS_2013_2014_2015.txt'
        # file_ecmwf = '/Volumes/Segon_HD/ecmwf/LaPalma/data_2013-2014-2015.txt'
        # tag_name = 'ecmwf_2013-2014-2015'  #'ecmwf_SOUTH_2013-2014-2015'
        # year_text = '2013-2014-2015'
        # epoch_text = 'winter'
        # label_ecmwf = 'ECMWF La Palma'
        # label_gdas  = 'GDAS La Pamla'
        # =========================


        self.Ns = Ns  # [cm^-3] molecular number density for standard air conditions
        self.ps = ps  # [hPa]   standard air pressure
        self.Ts = Ts  # [K]     standard air temperature
        self.Hs = Hs  # [km]    density scale hight for La Palma Winter
        # ----------------------------------------------------------------------

        # Plevel = np.array([1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700, 650, 600, 550, 500, 450, 400,
        #                    350, 300, 250, 225, 200, 175, 150, 125, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1])

        # MAGIC Winter parameters
        self.heightmw = heightmw
        self.n_mw = rhomw * 816.347
        self.func_magic = interp1d(self.heightmw * 1000., self.n_mw * np.exp(self.heightmw * 1000. / self.Hs),
                                   kind='cubic')

        # Prod3 Simulations (based on NRLMSISE)
        if self.observatory == 'north':
            prod3 = open(
                '/Volumes/Segon_HD/molecularprofiles/ecmwf/Prod3b_simulations_molecular_profiles/atmprof36_lapalma.dat')
            self.hprod3, nprod3, thickprod3, n1prod3, tempprod3, pprod3 = np.loadtxt(prod3, usecols=(0, 1, 2, 3, 4, 5),
                                                                                     unpack=True)
        elif self.observatory == 'shouth':
            prod3 = open(
                '/Volumes/Segon_HD/molecularprofiles/ecmwf/Prod3b_simulations_molecular_profiles/atmprof26_paranal.dat')
            self.hprod3, nprod3, thickprod3, n1prod3, tempprod3, pprod3 = np.loadtxt(prod3, usecols=(0, 1, 2, 3, 4, 5),
                                                                                     unpack=True)
        else:
            print('WRONG observatory. It must be "north" or "south" ')
            raise SystemExit

        self.n_prod3 = nprod3 * 816.347
        self.func_prod3 = interp1d(self.hprod3 * 1000., self.n_prod3 * np.exp(self.hprod3 * 1000. / self.Hs),
                                   kind='cubic')

    def retrievedata(self, model):
        """
        Function that reads one txt file containing the ecmwf or the gdas data and returns the data ready to produce
        the quality check plots

        :param model: valid options are: "ecmwf" or "gdas"
        :return: mjd_at_15km, month_at_15km, interpolated_density, int_density_at_15km, diff_with_magic, diff_with_prod3
        """
        print('This function is still not ready. Do not use it. Use instead retrievedatafromfiles')
        exit()

        # if model == 'ecmwf':
        #     mjd, year, month, day, hour, h, p, n = read_ecmwf(self.file_ecmwf, self.epoch_text)
        #     n = n / self.Ns
        #     mjd0, mjd_max = np.min(mjd), np.max(mjd)
        # if model == 'gdas':
        #     mjd, year, month, day, hour, h, n, p = read_gdas(self.file_gdas, self.epoch_text)
        #     year_ecmwf, month_ecmwf, day_ecmwf, hour_ecmwf = np.loadtxt(open(self.file_ecmwf), usecols=[1,2,3,4],
        #                                                                 unpack=True)
        #     mjd0    = date2mjd(year_ecmwf[0], month_ecmwf[0], day_ecmwf[0], hour_ecmwf[0])
        #     mjd_max = date2mjd(year_ecmwf[-1], month_ecmwf[-1], day_ecmwf[-1], hour_ecmwf[-1])
        #
        # self.x = np.linspace(1000., 25000., num=15, endpoint=True)
        # interpolated_density = []
        # diff_with_magic = []
        # diff_with_prod3 = []
        # mjd_i = mjd0
        # int_density_at_15km = []
        # mjd_at_15km = []
        # month_at_15km = []
        # steps = (mjd_max - mjd_i) / 0.25
        # print("Computing the extrapolation of the values of density for %s:" % (model) )
        # while mjd_i < (mjd_max + 0.25):
        #     sys.stdout.write('\r')
        #     k = int((mjd_i- mjd0)/0.25 * 100/steps)
        #     sys.stdout.write("[%-100s] %d%%" % ('=' * k, k))
        #     sys.stdout.write("")
        #     sys.stdout.flush()
        #     if np.all(mjd != mjd_i) or np.all(h[mjd == mjd_i] == 0.):
        #         mjd_i += 0.25
        #         continue
        #     func = interp1d(h[mjd == mjd_i], n[mjd == mjd_i] * np.exp(h[mjd == mjd_i] / self.Hs), kind='cubic')
        #     if np.isnan(func(self.x)).any():
        #         mjd_i += 0.25
        #         continue
        #
        #     interpolated_density.append(func(self.x))
        #     int_density_at_15km.append(func(self.x[8]))
        #     mjd_at_15km.append(mjd_i)
        #     month_at_15km.append(month[mjd == mjd_i])
        #     diff_with_magic.append((func(self.x) - self.func_magic(self.x)) / func(self.x))
        #     diff_with_prod3.append((func(self.x) - self.func_prod3(self.x)) / func(self.x))
        #     mjd_i += 0.25
        #
        # interpolated_density = np.asarray(interpolated_density)
        # int_density_at_15km = np.asarray(int_density_at_15km)
        # mjd_at_15km = np.asarray(mjd_at_15km)
        # month_at_15km = np.asarray(month_at_15km)
        # diff_with_magic = np.asarray(diff_with_magic)
        # diff_with_prod3 = np.asarray(diff_with_prod3)
        # print('\n')
        # return mjd_at_15km, month_at_15km, interpolated_density, int_density_at_15km, diff_with_magic, diff_with_prod3

    def retrieve_gdas_from_file(self):

        """
        Function that reads GDAS txt input files and returns quantities ready to plot
        :return: ()
        """

        self.mjd_gdas, self.year_gdas, self.month_gdas, self.day_gdas, self.hour_gdas, self.h_gdas, self.nns_gdas, \
        self.p_gdas = read_gdas(self.file_gdas, self.epoch_text)
        interpolated_density_gdas = []
        diff_gdas_with_magic = []
        diff_gdas_with_prod3 = []
        gdas_density_at_15km = []
        mjd_at_15km = []
        month_at_15km = []

        mjd = self.mjd_gdas[0]
        self.x = np.linspace(1000., 25000., num=15, endpoint=True)
        steps = (np.max(self.mjd_gdas) - mjd) / 0.25

        print("Computing the extrapolation of the values of density for GDAS:")
        while mjd < (np.max(self.mjd_gdas) + 0.25):
            # Percentage counter bar
            sys.stdout.write('\r')
            k = int((mjd - self.mjd_gdas[0]) / 0.25 * 100 / steps)
            sys.stdout.write("[%-100s] %d%%" % ('=' * k, k))
            sys.stdout.flush()
            # ---------------------------
            if np.all(self.mjd_gdas != mjd) or np.all(self.h_gdas[self.mjd_gdas == mjd] == 0.):
                mjd += 0.25
                continue
            func_gdas = interp1d(self.h_gdas[self.mjd_gdas == mjd], self.nns_gdas[self.mjd_gdas == mjd] *
                                 np.exp(self.h_gdas[self.mjd_gdas == mjd] / self.Hs), kind='cubic')
            if np.isnan(func_gdas(self.x)).any():
                mjd += 0.25
                continue

            interpolated_density_gdas.append(func_gdas(self.x))
            gdas_density_at_15km.append(func_gdas(self.x[8]))
            mjd_at_15km.append(mjd)
            month_at_15km.append(self.month_gdas[self.mjd_gdas == mjd])
            diff_gdas_with_magic.append((func_gdas(self.x) - self.func_magic(self.x)) / func_gdas(self.x))
            diff_gdas_with_prod3.append((func_gdas(self.x) - self.func_prod3(self.x)) / func_gdas(self.x))
            mjd += 0.25

        print('\n')

        self.interpolated_density_gdas = np.asarray(interpolated_density_gdas)
        self.gdas_density_at_15km = np.asarray(gdas_density_at_15km)
        self.diff_gdas_with_magic = np.asarray(diff_gdas_with_magic)
        self.diff_gdas_with_prod3 = np.asarray(diff_gdas_with_prod3)
        self.mjd_at_15km_gdas = np.asarray(mjd_at_15km)
        self.gdas_averages = compute_averages_std(self.interpolated_density_gdas)
        # Differences w.r.t. MAGIC W model
        self.diff_gdas_MAGIC = compute_averages_std(self.diff_gdas_with_magic)
        # Differences w.r.t. Prod3
        self.diff_gdas_PROD3 = compute_averages_std(self.diff_gdas_with_prod3)

    def retrieve_ecmwf_from_file(self):

        """
        Function that reads ECMWF txt input files and returns quantities ready to plot
        :return: ()
        """

        if not os.path.exists((self.file_ecmwf).split('.')[0] + '.txt'):
            grib_file = self.file_ecmwf
            readgribfile2text(grib_file, self.year, self.observatory)

        self.mjd_ecmwf, self.year_ecmwf, self.month_ecmwf, self.day_ecmwf, self.hour_ecmwf, self.h_ecmwf, self.p_ecmwf, \
        self.n_ecmwf = read_ecmwf(self.file_ecmwf.split('.')[0] + '.txt', self.epoch_text)

        interpolated_density_ecmwf = []
        diff_ecmwf_with_magic = []
        diff_ecmwf_with_prod3 = []
        ecmwf_density_at_15km = []
        mjd_at_15km = []
        month_at_15km = []

        mjd = self.mjd_ecmwf[0]
        self.x = np.linspace(1000., 25000., num=15, endpoint=True)
        steps = (np.max(self.mjd_ecmwf) - mjd) / 0.25

        print("Computing the extrapolation of the values of density for ECMWF:")
        while mjd < (np.max(self.mjd_ecmwf) + 0.25):
            # Percentage counter bar
            sys.stdout.write('\r')
            k = int((mjd - self.mjd_ecmwf[0]) / 0.25 * 100 / steps)
            sys.stdout.write("[%-100s] %d%%" % ('=' * k, k))
            sys.stdout.flush()
            # ---------------------------
            func_ecmwf = interp1d(self.h_ecmwf[self.mjd_ecmwf == mjd], self.n_ecmwf[self.mjd_ecmwf == mjd] / self.Ns *
                                  np.exp(self.h_ecmwf[self.mjd_ecmwf == mjd] / self.Hs), kind='cubic')

            interpolated_density_ecmwf.append(func_ecmwf(self.x))
            ecmwf_density_at_15km.append(func_ecmwf(self.x[8]))
            mjd_at_15km.append(mjd)
            month_at_15km.append(self.month_ecmwf[self.mjd_ecmwf == mjd])
            diff_ecmwf_with_magic.append((func_ecmwf(self.x) - self.func_magic(self.x)) / func_ecmwf(self.x))
            diff_ecmwf_with_prod3.append((func_ecmwf(self.x) - self.func_prod3(self.x)) / func_ecmwf(self.x))

            mjd += 0.25

        print('\n')
        self.interpolated_density_ecmwf = np.asarray(interpolated_density_ecmwf)
        self.ecmwf_density_at_15km = np.asarray(ecmwf_density_at_15km)
        self.mjd_at_15km_ecmwf = np.asarray(mjd_at_15km)
        self.diff_ecmwf_with_magic = np.asarray(diff_ecmwf_with_magic)
        self.diff_ecmwf_with_prod3 = np.asarray(diff_ecmwf_with_prod3)
        self.ecmwf_averages = compute_averages_std(self.interpolated_density_ecmwf)
        # Differences w.r.t. MAGIC W model
        self.diff_ecmwf_MAGIC = compute_averages_std(self.diff_ecmwf_with_magic)
        # Differences w.r.t. Prod3
        self.diff_ecmwf_PROD3 = compute_averages_std(self.diff_ecmwf_with_prod3)

        # def OLD_retrieve_data_from_file(self):
        #
        #     """
        #     Function that reads both ECMWF and GDAS txt input files and returns quantities ready to plot
        #     :return: ()
        #     """
        #
        #     if not os.path.exists(self.file_ecmwf):
        #         readgribfile2text(self.grib_file, self.year, self.observatory)
        #
        #     self.mjd_ecmwf, self.year_ecmwf, self.month_ecmwf, self.day_ecmwf, self.hour_ecmwf, self.h_ecmwf, self.p_ecmwf, \
        #     self.n_ecmwf = read_ecmwf(self.file_ecmwf, self.epoch_text)
        #
        #     interpolated_density_ecmwf = []
        #     interpolated_density_gdas = []
        #     diff_ecmwf_with_magic = []
        #     diff_ecmwf_with_prod3 = []
        #     mjd = self.mjd_ecmwf[0]
        #     ecmwf_density_at_15km = []
        #     mjd_at_15km = []
        #     month_at_15km = []
        #
        #     if self.file_gdas:
        #         self.mjd_gdas, self.year_gdas, self.month_gdas, self.day_gdas, self.hour_gdas, self.h_gdas, self.nns_gdas, \
        #         self.p_gdas = read_gdas(self.file_gdas, self.epoch_text)
        #         interpolated_density_gdas = []
        #         diff_gdas_with_magic = []
        #         diff_gdas_with_prod3 = []
        #         gdas_density_at_15km = []
        #
        #     self.x = np.linspace(1000., 25000., num=15, endpoint=True)
        #     steps = (np.max(self.mjd_ecmwf) - mjd) / 0.25
        #     print("Computing the extrapolation of the values of density for GDAS and ECMWF:")
        #     while mjd < (np.max(self.mjd_ecmwf) + 0.25):
        #         # Percentage counter bar
        #         sys.stdout.write('\r')
        #         k = int((mjd - self.mjd_ecmwf[0]) / 0.25 * 100 / steps)
        #         sys.stdout.write("[%-100s] %d%%" % ('=' * k, k))
        #         sys.stdout.flush()
        #         # ---------------------------
        #         if self.file_gdas:
        #             if np.all(self.mjd_gdas != mjd) or np.all(self.h_gdas[self.mjd_gdas == mjd] == 0.):
        #                 mjd += 0.25
        #                 continue
        #             func_gdas = interp1d(self.h_gdas[self.mjd_gdas == mjd], self.nns_gdas[self.mjd_gdas == mjd] *
        #                                  np.exp(self.h_gdas[self.mjd_gdas == mjd] / self.Hs), kind='cubic')
        #             if np.isnan(func_gdas(self.x)).any():
        #                 mjd += 0.25
        #                 continue
        #
        #             interpolated_density_gdas.append(func_gdas(self.x))
        #             gdas_density_at_15km.append(func_gdas(self.x[8]))
        #             diff_gdas_with_magic.append((func_gdas(self.x) - self.func_magic(self.x)) / func_gdas(self.x))
        #             diff_gdas_with_prod3.append((func_gdas(self.x) - self.func_prod3(self.x)) / func_gdas(self.x))
        #
        #         func_ecmwf = interp1d(self.h_ecmwf[self.mjd_ecmwf == mjd],
        #                               self.n_ecmwf[self.mjd_ecmwf == mjd] / self.Ns *
        #                               np.exp(self.h_ecmwf[self.mjd_ecmwf == mjd] / self.Hs), kind='cubic')
        #
        #         interpolated_density_ecmwf.append(func_ecmwf(self.x))
        #         ecmwf_density_at_15km.append(func_ecmwf(self.x[8]))
        #         mjd_at_15km.append(mjd)
        #         month_at_15km.append(self.month_ecmwf[self.mjd_ecmwf == mjd])
        #         diff_ecmwf_with_magic.append((func_ecmwf(self.x) - self.func_magic(self.x)) / func_ecmwf(self.x))
        #         diff_ecmwf_with_prod3.append((func_ecmwf(self.x) - self.func_prod3(self.x)) / func_ecmwf(self.x))
        #
        #         mjd += 0.25
        #
        #     print('\n')
        #     self.interpolated_density_ecmwf = np.asarray(interpolated_density_ecmwf)
        #     self.ecmwf_density_at_15km = np.asarray(ecmwf_density_at_15km)
        #     self.mjd_at_15km_ecmwf = np.asarray(mjd_at_15km)
        #     self.diff_ecmwf_with_magic = np.asarray(diff_ecmwf_with_magic)
        #     self.diff_ecmwf_with_prod3 = np.asarray(diff_ecmwf_with_prod3)
        #
        #     # self.mjd_at_15km_gdas, self.month_at_15km_gdas, self.interpolated_density_gdas, self.gdas_density_at_15km, \
        #     # self.diff_gdas_with_magic, self.diff_gdas_with_prod3 = self.retrievedata('gdas')
        #     #
        #     # self.mjd_at_15km_ecmwf, self.month_at_15km_ecmwf, self.interpolated_density_ecmwf, self.ecmwf_density_at_15km, \
        #     # self.diff_ecmwf_with_magic, self.diff_ecmwf_with_prod3 = self.retrievedata('ecmwf')
        #
        #     if self.file_gdas:
        #         self.interpolated_density_gdas = np.asarray(interpolated_density_gdas)
        #         self.gdas_density_at_15km = np.asarray(gdas_density_at_15km)
        #         self.diff_gdas_with_magic = np.asarray(diff_gdas_with_magic)
        #         self.diff_gdas_with_prod3 = np.asarray(diff_gdas_with_prod3)
        #         print('relative difference between ECMWF and GDAS means= %5.2f%%' % (
        #         (np.mean(self.ecmwf_density_at_15km) -
        #          np.mean(self.gdas_density_at_15km)) /
        #         np.mean(self.ecmwf_density_at_15km) * 100.))
        #         self.gdas_averages = compute_averages_std(self.interpolated_density_gdas)
        #         self.differences_avgs = compute_averages_std((self.interpolated_density_gdas -
        #                                                       self.interpolated_density_ecmwf) /
        #                                                      self.interpolated_density_gdas)
        #         # Differences w.r.t. MAGIC W model
        #         self.diff_gdas_MAGIC = compute_averages_std(self.diff_gdas_with_magic)
        #         # Differences w.r.t. Prod3
        #         self.diff_gdas_PROD3 = compute_averages_std(self.diff_gdas_with_prod3)
        #
        #     self.ecmwf_averages = compute_averages_std(self.interpolated_density_ecmwf)
        #     # self.differences_avgs = compute_averages_std_simple((self.ecmwf_averages[0] - self.gdas_averages[0])
        #     #                                             / self.ecmwf_averages[0])
        #     # Differences w.r.t. MAGIC W model
        #     self.diff_ecmwf_MAGIC = compute_averages_std(self.diff_ecmwf_with_magic)
        #     # Differences w.r.t. Prod3
        #     self.diff_ecmwf_PROD3 = compute_averages_std(self.diff_ecmwf_with_prod3)

    def difference_bteween_gdas_ecmwf_(self):
        print('relative difference between ECMWF and GDAS means= %5.2f%%' % ((np.mean(self.ecmwf_density_at_15km) -
                                                                              np.mean(self.gdas_density_at_15km)) /
                                                                             np.mean(
                                                                                 self.ecmwf_density_at_15km) * 100.))
        if len(self.interpolated_density_ecmwf) > len(self.interpolated_density_gdas):
            self.interpolated_density_ecmwf = self.interpolated_density_ecmwf[:len(self.interpolated_density_gdas)]
        else:
            self.interpolated_density_gdas = self.interpolated_density_gdas[:len(self.interpolated_density_ecmwf)]

            self.differences_avgs = compute_averages_std((self.interpolated_density_gdas -
                                                          self.interpolated_density_ecmwf) /
                                                         self.interpolated_density_gdas)

    # =======================================================================================================
    def plot_average_at_15km(self):
        """
        Function that produces a plot of the averaged density at 15 km for both GDAS and ECMWF data
        :return:
        """

        # Plotting the differences between the averages of GDAS and ECMWF data for each day, in bins of height.
        print('plotting data at 15 km in 1-day bins:')

        fig = plt.figure()
        ax = fig.add_subplot(111)
        try:
            self.gdas_density_at_15km
            ax.plot(self.mjd_at_15km_gdas, self.gdas_density_at_15km, 'o', color='#FF9999', markersize=0.75,
                    label=self.label_gdas + ' ' + self.observatory, alpha=0.3)
            ax.hlines(np.mean(self.gdas_density_at_15km), np.min(self.mjd_at_15km_ecmwf) - 10.,
                      np.max(self.mjd_at_15km_ecmwf) + 10.,
                      colors='r', linestyle='solid', lw=2., zorder=10)
        except:
            print('No GDAS file loaded')

        ax.plot(self.mjd_at_15km_ecmwf, self.ecmwf_density_at_15km, 'o', color='#99CCFF', markersize=0.75,
                label=self.label_ecmwf + ' ' + self.observatory, alpha=0.3)

        # Make a function to plot vertical lines on starting of different years and pass it as option to the main function call
        # o = 0
        # years = self.year
        # for y in [2013, 2014, 2015]:
        #     ax.vlines(np.min(self.mjd_at_15km_ecmwf) + 365.25 * o, 0.775, 1.0, color='0.7', linewidth=1.)
        #     ax.text(np.min(self.mjd_at_15km_ecmwf) + 365.25 * o - 7.25, 0.77, str(y), rotation=90, color='0.7')
        #     ax.vlines(np.min(self.mjd_at_15km_ecmwf) + 365.25 * o + 181., 0.6, 1.0, color='0.7', linestyles='dotted',
        #               linewidth=1.)
        #     o += 1

        ax.hlines(np.mean(self.ecmwf_density_at_15km), np.min(self.mjd_at_15km_ecmwf) - 10.,
                  np.max(self.mjd_at_15km_ecmwf) + 10.,
                  colors='b', linestyle='solid', lw=2., zorder=10)

        ax.legend(loc='best', numpoints=1)
        ax.set_xlabel('MJD')
        ax.set_ylabel('$n_{\\rm day}/N_{\\rm s} \\cdot e^{(h/H_{\\rm s})}$')
        ax.set_ylim(0.757, 0.939)
        ax.set_title('Density over time at h = 15 km (' + str(self.epoch_text) + ')')
        fig.savefig(self.output_plot_name + '_at_15_km.eps', bbox_inches='tight')
        fig.savefig(self.output_plot_name + '_at_15_km.png', bbox_inches='tight', dpi=300)

    # fig.clf()

    # gdas_density_at_15km  = gdas_density_at_15km.reshape(-1,7)
    # ecmwf_density_at_15km = ecmwf_density_at_15km.reshape(-1,7)
    # mjd_at_15km_ecmwf   = mjd_at_15km_ecmwf.reshape(-1,7)
    # avg_gdas_at_15    = []
    # std_gdas_at_15    = []
    # avg_ecmwf_at_15   = []
    # std_ecmwf_at_15   = []
    # mjd_at_15         = []

    # for i in np.arange(len(gdas_density_at_15km[0,:])):
    #   avg_gdas_at_15.append(np.average(gdas_density_at_15km[i,:]))
    #   avg_ecmwf_at_15.append(np.average(ecmwf_density_at_15km[i,:]))
    #   std_gdas_at_15.append(robust.mad(gdas_density_at_15km[i,:]))
    #   std_ecmwf_at_15.append(robust.mad(ecmwf_density_at_15km[i,:]))

    # fig = plt.figure()
    # ax  = fig.add_subplot(111)
    # ax.plot(mjd_at_15km_ecmwf, gdas_density_at_15km, 'or', markersize=1.5, label='GDAS')
    # ax.plot(mjd_at_15km_ecmwf, ecmwf_density_at_15km, 'sb', markersize=1.5, label='ECMWF')
    # ax.legend(loc='best', numpoints=1)
    # ax.set_xlabel('MJD')
    # ax.set_ylabel('$n_{\\rm day}/N_{\\rm s} \\cdot e^{(h/H_{\\rm s})}$')
    # fig.savefig(tag_name+'_at_15_km.eps', bbox_inches='tight')
    # fig.savefig(tag_name+'_at_15_km.png', bbox_inches='tight', dpi=300)
    # fig.clf()

    def plot_differences_averages(self):
        """
        Function that provides a plot of the relative difference between the ECMWF and GDAS density models
        :return:
        """
        try:
            self.differences_avgs
            print('Computing the averages, std dev and peak-to-peak values for the differences:')
            print('plotting averaged data values for selected epoch')
            fig = plt.figure()
            ax = fig.add_subplot(111)

            eb1 = ax.errorbar(self.x, self.differences_avgs[0],
                              yerr=[self.differences_avgs[3], self.differences_avgs[2]],
                              fmt='o', color='b', capsize=0.5, mec='b', ms=1.)
            eb1[-1][0].set_linestyle(':')
            ax.errorbar(self.x, self.differences_avgs[0], yerr=self.differences_avgs[1], fmt=':', color='g',
                        elinewidth=2.5)

            ax.set_title('Relative Difference GDAS-ECMWF %s %s' % (self.year, self.epoch_text))
            ax.set_xlabel('h a.s.l. [m]')
            ax.set_ylabel('Rel. Difference (GDAS - ECMWF)')
            ax.set_xlim(0., 25100.)
            ax.set_ylim(-0.11, 0.09)
            ax.xaxis.set_minor_locator(MultipleLocator(1000))
            ax.xaxis.set_major_locator(MultipleLocator(2000))
            ax.yaxis.set_major_locator(MultipleLocator(0.02))
            # ax.legend(loc='best', numpoints=1)
            ax.grid(which='both', axis='y', color='0.8')
            fig.savefig('differences_' + self.output_plot_name + '.eps', bbox_inches='tight')
            fig.savefig('differences_' + self.output_plot_name + '.png', bbox_inches='tight', dpi=300)
        except:
            print('Cannot plot differences and averages now...')

    def plot_differences_wrt_magic(self):
        """
        Function that plots the difference between the ECMWF and GDAS density models w.r.t the MAGIC Winter model
        :return:
        """

        print('Computing the averages, std dev and peak-to-peak values for the differences wrt MAGIC Winter model:')
        print('plotting averaged data values for selected epoch')
        fig = plt.figure()
        ax = fig.add_subplot(111)

        try:
            self.diff_gdas_MAGIC[0]
            eb1 = ax.errorbar(self.x, self.diff_gdas_MAGIC[0], yerr=[self.diff_gdas_MAGIC[3], self.diff_gdas_MAGIC[2]],
                              fmt='o', color='r', capsize=0.5, mec='r', ms=1., label=self.label_gdas)
            eb1[-1][0].set_linestyle(':')
            ax.errorbar(self.x, self.diff_gdas_MAGIC[0], yerr=self.diff_gdas_MAGIC[1], fmt=':', color='r',
                        elinewidth=3.1)
        except:
            print('No GDAS file loaded')

        eb2 = ax.errorbar(self.x + 175., self.diff_ecmwf_MAGIC[0], yerr=[self.diff_ecmwf_MAGIC[3],
                                                                         self.diff_ecmwf_MAGIC[2]], fmt='o', color='b',
                          capsize=0.5, mec='b', ms=1., label=self.label_ecmwf)
        eb2[-1][0].set_linestyle(':')
        ax.errorbar(self.x + 175., self.diff_ecmwf_MAGIC[0], yerr=self.diff_ecmwf_MAGIC[1], fmt=':', color='b',
                    elinewidth=3.1)

        ax.set_title('Relative Difference w.r.t MAGIC W model %s %s' % (self.year, self.epoch_text))
        ax.set_xlabel('h a.s.l. [m]')
        ax.set_ylabel('Rel. Difference (model - MW)')
        ax.set_xlim(0., 25100.)
        ax.set_ylim(-0.11, 0.09)
        ax.xaxis.set_minor_locator(MultipleLocator(1000))
        ax.xaxis.set_major_locator(MultipleLocator(2000))
        ax.yaxis.set_major_locator(MultipleLocator(0.02))
        ax.legend(loc='best', numpoints=1)
        ax.grid(which='both', axis='y', color='0.8')
        fig.savefig('differences_wrt_MAGIC_' + self.output_plot_name + '.eps', bbox_inches='tight')
        fig.savefig('differences_wrt_MAGIC_' + self.output_plot_name + '.png', bbox_inches='tight', dpi=300)

    def plot_differences_wrt_other(self, model):
        print('Computing the averages, std dev and peak-to-peak values for the differences wrt MAGIC Winter model:')
        print('plotting averaged data values for selected epoch')
        print('NEED TO FINisH THIS!')
        if model == 'MW':
            try:
                self.diff_gdas_MAGIC[0]
                diff_gdas = self.diff_gdas_MAGIC[0]
                ediff_gdas_pp = [self.diff_gdas_MAGIC[3], self.diff_gdas_MAGIC[2]]
                ediff_gdas = self.diff_gdas_MAGIC[1]
            except:
                print('No GDAS file loaded')
            diff_ecmwf = self.diff_ecmwf_MAGIC[0]
            ediff_ecmwf_pp = [self.diff_ecmwf_MAGIC[3], self.diff_ecmwf_MAGIC[2]]
            ediff_ecmwf = self.diff_ecmwf_MAGIC[1]
        elif model == 'PROD3':
            if self.diff_gdas_PROD3[0].all:
                diff_gdas = self.diff_gdas_PROD3[0]
                ediff_gdas_pp = [self.diff_gdas_PROD3[3], self.diff_gdas_PROD3[2]]
                ediff_gdas = self.diff_gdas_MAGIC[1]
            diff_ecmwf = self.diff_ecmwf_PROD3[0]
            ediff_ecmwf_pp = [self.diff_ecmwf_PROD3[3], self.diff_ecmwf_PROD3[2]]
            ediff_ecmwf = self.diff_ecmwf_PROD3[1]
        else:
            return (print('Wrong model name'))

        fig = plt.figure()
        ax = fig.add_subplot(111)

        try:
            diff_gdas
            eb1 = ax.errorbar(self.x, diff_gdas, yerr=ediff_gdas_pp, fmt='o', color='r', capsize=0.5, mec='r', ms=1.,
                              label=self.label_gdas)
            eb1[-1][0].set_linestyle(':')
            ax.errorbar(self.x, diff_gdas, yerr=ediff_gdas, fmt=':', color='r', elinewidth=3.1)
        except:
            print('No GDAS file loaded')

        eb2 = ax.errorbar(self.x + 175., diff_ecmwf, yerr=ediff_ecmwf_pp, fmt='o', color='b', capsize=0.5, mec='b',
                          ms=1., label=self.label_ecmwf)
        eb2[-1][0].set_linestyle(':')
        ax.errorbar(self.x + 175., diff_ecmwf, yerr=ediff_ecmwf, fmt=':', color='b',
                    elinewidth=3.1)

        ax.set_title('Relative Difference w.r.t %s model, year: %s, epoch: %s' % (model, self.year,
                                                                                  self.epoch_text))
        ax.set_xlabel('h a.s.l. [m]')
        ax.set_ylabel('Rel. Difference')
        ax.set_xlim(0., 25100.)
        ax.set_ylim(-0.11, 0.09)
        ax.xaxis.set_minor_locator(MultipleLocator(1000))
        ax.xaxis.set_major_locator(MultipleLocator(2000))
        ax.yaxis.set_major_locator(MultipleLocator(0.02))
        ax.legend(loc='best', numpoints=1)
        ax.grid(which='both', axis='y', color='0.8')
        fig.savefig('differences_wrt_PROD3_' + self.output_plot_name + '.eps', bbox_inches='tight')
        fig.savefig('differences_wrt_PROD3_' + self.output_plot_name + '.png', bbox_inches='tight', dpi=300)

    def plot_models_comparison(self, model):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        try:
            self.gdas_averages[0]
            eb1 = ax.errorbar(self.x, self.gdas_averages[0], yerr=[self.gdas_averages[3], self.gdas_averages[2]],
                              fmt='o',
                              color='r', capsize=0.5, mec='r', ms=1., label=self.label_gdas)
            eb1[-1][0].set_linestyle(':')
            ax.errorbar(self.x, self.gdas_averages[0], yerr=self.gdas_averages[1], fmt=':', color='r', elinewidth=3.)
        except:
            print('No GDAS file loaded')

        eb2 = ax.errorbar(self.x + 175., self.ecmwf_averages[0], yerr=[self.ecmwf_averages[3], self.ecmwf_averages[2]],
                          fmt='o', color='b', capsize=0.5, mec='b', ms=1., label=self.label_ecmwf)
        eb2[-1][0].set_linestyle(':')
        ax.errorbar(self.x + 175., self.ecmwf_averages[0], yerr=self.ecmwf_averages[1], fmt=':', color='b',
                    elinewidth=3.)

        if model == 'MW':
            ax.plot(self.heightmw * 1000., self.n_mw * np.exp(self.heightmw * 1000. / self.Hs), '-', color='grey',
                    label='MAGIC W')
        elif model == 'PROD3':
            ax.plot(self.hprod3 * 1000., self.n_prod3 * np.exp(self.hprod3 * 1000. / self.Hs), '-', color='0.2',
                    label='Prod3 ' + self.observatory)
        elif model == 'both':
            ax.plot(self.heightmw * 1000., self.n_mw * np.exp(self.heightmw * 1000. / self.Hs), '-', color='grey',
                    label='MAGIC W')
            ax.plot(self.hprod3 * 1000., self.n_prod3 * np.exp(self.hprod3 * 1000. / self.Hs), '-', color='0.2',
                    label='Prod3 ' + self.observatory)

        ax.set_title(self.label_ecmwf + ' ' + str(self.year) + ' ' + str(self.epoch_text))
        ax.set_xlabel('h a.s.l. [m]')
        ax.set_ylabel('$n_{\\rm day}/N_{\\rm s} \\cdot e^{(h/H_{\\rm s})}$')
        ax.set_xlim(0., 25100.)
        ax.set_ylim(0.4, 1.2)
        ax.xaxis.set_minor_locator(MultipleLocator(1000))
        ax.xaxis.set_major_locator(MultipleLocator(2000))
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
        ax.legend(loc='best', numpoints=1)
        ax.grid(which='both', axis='y', color='0.8')
        fig.savefig('comparison_' + self.output_plot_name + '.eps', bbox_inches='tight')
        fig.savefig('comparison_' + self.output_plot_name + '.png', bbox_inches='tight', dpi=300)

    def print_to_text_file(self):
        textfile = open(self.output_plot_name + '_to_text_file.txt', 'w')

        try:
            self.gdas_averages[0]
            print('# x[i], gdas_avg, gdas_rms, gdas_p2p_m, gdas_p2p_p, ecmwf_avg, ecmwf_rms, ecmwf_p2p_m, ecmwf_p2p_p',
                  file=textfile)
            for i in np.arange(len(self.x)):
                print(self.x[i], self.gdas_averages[0][i], self.gdas_averages[1][i], self.gdas_averages[3][i],
                      self.gdas_averages[2][i], self.ecmwf_averages[0][i], self.ecmwf_averages[1][i],
                      self.ecmwf_averages[3][i],
                      self.ecmwf_averages[2][i], self.differences_avgs[0][i], self.differences_avgs[1][i],
                      self.differences_avgs[3][i], self.differences_avgs[2][i], self.diff_gdas_MAGIC[0][i],
                      self.diff_gdas_MAGIC[1][i], self.diff_gdas_MAGIC[3][i], self.diff_gdas_MAGIC[2][i],
                      self.diff_ecmwf_MAGIC[0][i], self.diff_ecmwf_MAGIC[1][i], self.diff_ecmwf_MAGIC[3][i],
                      self.diff_ecmwf_MAGIC[2][i], file=textfile)
            textfile.close()
        except:
            print('No GDAS file loaded')

            print('# x[i], ecmwf_avg, ecmwf_rms, ecmwf_p2p_m, ecmwf_p2p_p',
                  file=textfile)
            for i in np.arange(len(self.x)):
                print(self.x[i], self.ecmwf_averages[0][i], self.ecmwf_averages[1][i], self.ecmwf_averages[3][i],
                      self.ecmwf_averages[2][i], self.diff_ecmwf_MAGIC[0][i], self.diff_ecmwf_MAGIC[1][i],
                      self.diff_ecmwf_MAGIC[3][i], self.diff_ecmwf_MAGIC[2][i], file=textfile)
            textfile.close()

