import os
import os.path
import sys
import matplotlib.pyplot as plt
from tqdm import *
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import interp1d
from molecularprofiles.utils.grib_utils import *
from molecularprofiles.utils.read_txt_gribfile import *
from molecularprofiles.utils.plot_settings import settings
from molecularprofiles.utils.magic_winter_profile import heightmw, rhomw
from molecularprofiles.utils.meteorological_constants import *
from molecularprofiles.utils.humidity import *

settings()


class MolecularProfile:
    def __init__(self, data_file, data_server='ECMWF', tag_name='myplots', observatory='north'):

        """
        This class provides with a series of functions to analyze the quality of the data for both
        CTA-North and CTA-South.

        :param data_file: txt file containing the data (string)
        :param tag_name: name to be given to the output files (string)
        :param epoch_text: season to be analyzed. Valid values are: "all", "winter", "summer", "intermediate" (string) 
        :param data_server: label to be put in some of the output plots (string)
        :param observatory: valid options are: "north", "south"

        Functions within this class:

        retrievedatafromfile: same as retrievedata but in this case it also performs the calculation of the averages
                              for several variables
        plot_average_at_15km
        plot_differences_wrt_magic
        plot_differences_wrt_prod3
        plot_models_comparison
        print_to_text_file

        """

        self.data_file = data_file
        self.tag_name = tag_name
        self.data_server = data_server
        self.observatory = observatory

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

    def _get_prod3sim_data(self):

        # TODO: change the directory definition. Make an import since everything is now recognised as package

        MOLECULARPROFILES_DIR = os.environ.get('MOLECULARPROFILES_DIR')
        PROD3_DIR = MOLECULARPROFILES_DIR + 'molecularprofiles/Prod3b_simulations/'
        # Prod3 Simulations (based on NRLMSISE)
        if self.observatory == 'north':
            prod3 = open(PROD3_DIR + '/atmprof36_lapalma.dat')
            self.hprod3, nprod3, thickprod3, n1prod3, tempprod3, pprod3 = np.loadtxt(prod3, usecols=(0, 1, 2, 3, 4, 5),
                                                                                     unpack=True)
        elif self.observatory == 'shouth':
            prod3 = open(PROD3_DIR + '/atmprof26_paranal.dat')
            self.hprod3, nprod3, thickprod3, n1prod3, tempprod3, pprod3 = np.loadtxt(prod3, usecols=(0, 1, 2, 3, 4, 5),
                                                                                     unpack=True)
        else:
            print('WRONG observatory. It must be "north" or "south" ')
            raise SystemExit

        self.n_prod3 = nprod3 * 816.347
        self.func_prod3 = interp1d(self.hprod3 * 1000., self.n_prod3 * np.exp(self.hprod3 * 1000. / self.Hs),
                                   kind='cubic')

    def get_data(self, epoch='all'):

        """
        Function that reads ECMWF or GDAS txt input files and returns quantities ready to plot.
        If the input filename does not exist, the program searches for the same
        input name with .grib extension and extracts the data from it

        Input: epoch_text: (str) can be "winter", "summer", "intermediate", "all". Default is "all"

        :return:
            self.date (YYYYMMDD)
            self.year (YYYY)
            self.month (MM)
            self.day (DD)
            self.hour (HH)
            self.mjd (DDDD.D)
            self.p (hPa)
            self.h (m)
            self.n (cm^-3)
            self.T (K)
            self.U (m s^-1)
            self.V (m s^-1)
            self.RH (%)
            self.interpolated_density (cm^-3)
            self.x (m)
            self.density_at_15km (cm^-3)
            self.mjd_at_15km (DDDD.D)
            self.averages (
        """

        if not os.path.exists((self.data_file).split('.')[0] + '.txt'):
            grib_file = self.data_file
            readgribfile2text(grib_file, self.observatory, gridstep=0.75)

        self.epoch_text = epoch
        self.output_plot_name = self.tag_name + '_' + self.epoch_text

        self.date, self.year, self.month, self.day, self.hour, self.mjd, self.p, self.h, \
        self.n, self.T, self.U, self.V, self.RH = read_file(self.data_file.split('.')[0] + '.txt', self.epoch_text)

        interpolated_density = []
        density_at_15km = []
        mjd_at_15km = []
        self.x = np.linspace(2200., 25000., num=15, endpoint=True)

        print("Computing the extrapolation of the values of density:")
        print("(This is to make it easier to compare ECMWF and GDAS, or any other")
        print("weather model)")
        pbar = tqdm(total=len(np.unique(self.mjd)))

        for mjd in np.unique(self.mjd):
            pbar.update(1)
            func = interp1d(self.h[self.mjd == mjd], self.n[self.mjd == mjd] / self.Ns
                                  * np.exp(self.h[self.mjd == mjd] / self.Hs), kind='cubic')

            int_dens = func(self.x)
            interpolated_density.append(int_dens)
            density_at_15km.append(func(self.x[8]))
            mjd_at_15km.append(mjd)
        pbar.close()

        print('\n')
        self.interpolated_density = np.asarray(interpolated_density)
        self.density_at_15km = np.asarray(density_at_15km)
        self.mjd_at_15km = np.asarray(mjd_at_15km)
        self.averages = compute_averages_std(self.interpolated_density)

    def compute_mass_density(self, air='moist'):
        """
        Uses the functions DensityMoistAir, MolarFractionWaterVapor and Compressibility from the LIDAR_analysis module humidity.py
        input:
            air: (str, optional) must be 'moist' or 'dry'

        :return: density [kg/m^3], std_dev(density) [kg/m^3], peak2peak_minus, peak2peak_plus
        """
        C = 415 # CO2 average global concentration in ppm
        Z = []
        rho = []
        Xw = []
        interpolated_rho = []

        pbar = tqdm(total=len(self.p))
        for i in np.arange(len(self.p)):
            pbar.update(1)
            if air == 'moist':
                Xw.append(MolarFractionWaterVapor(self.p[i], self.T[i], self.RH[i]))
                Z.append(Compressibility(self.p[i], self.T[i], Xw[i]))
                rho.append(DensityMoistAir(self.p[i] * 100., self.T[i], Z[i], Xw[i], C))

            elif air == 'dry':
                Z.append(Compressibility(self.p[i], self.T[i], 0.0))
                rho.append(DensityMoistAir(self.p[i]*100., self.T[i], Z[i], 0.0, C))
            else:
                print('Wrong air condition. It must be "moist" or "dry". Aborting!')
                sys.exit()
        pbar.close()

        rho = np.asarray(rho)
        pbar = tqdm(total=len(np.unique(self.mjd)))
        for mjd in np.unique(self.mjd):
            pbar.update(1)
            func_rho = interp1d(self.h[self.mjd == mjd], rho[self.mjd == mjd] *
                                np.exp(self.h[self.mjd == mjd] / self.Hs), kind = 'cubic', fill_value='extrapolate')
            interpolated_rho.append(func_rho(self.x))
        pbar.close()
        interpolated_rho = np.asarray(interpolated_rho)
        rho, e_rho, rho_pp, rho_pm = compute_averages_std(interpolated_rho)
        return rho, e_rho, rho_pp, rho_pm

    def compute_diff_wrt_model(self):

        diff_with_magic = []
        diff_with_prod3 = []

        self._get_prod3sim_data()
        self.x = np.linspace(1000., 25000., num=15, endpoint=True)

        print("Computing the differences of the values of density:")
        for i in tqdm(np.arange(len(self.interpolated_density))):
            # Percentage counter bar
            # sys.stdout.write('\r')
            # k = 100
            # sys.stdout.write("[%-100s] %d%%" % ('=' * k, k))
            # sys.stdout.flush()
            # ---------------------------
            diff_with_magic.append((self.interpolated_density[i] - self.func_magic(self.x))
                                         / self.interpolated_density[i])
            diff_with_prod3.append((self.interpolated_density[i] - self.func_prod3(self.x))
                                         / self.interpolated_density[i])

        self.diff_with_magic = np.asarray(diff_with_magic)
        self.diff_with_prod3 = np.asarray(diff_with_prod3)
        self.diff_MAGIC = compute_averages_std(self.diff_with_magic)
        self.diff_PROD3 = compute_averages_std(self.diff_with_prod3)

    # =======================================================================================================

    def plot_average_at_15km(self):
        """
        Function that produces a plot of the averaged density at 15 km
        :return: 
        """

        print('plotting data at 15 km in 1-day bins:')

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.mjd_at_15km, self.density_at_15km, 'o', color='#99CCFF', markersize=1.2,
                label=self.data_server + ' ' + self.observatory, alpha=0.8)

        ax.legend(loc='best', numpoints=1)
        ax.set_xlabel('MJD')
        ax.set_ylabel('$n_{\\rm day}/N_{\\rm s} \\cdot e^{(h/H_{\\rm s})}$')
        ax.set_ylim(np.min(self.density_at_15km) * 0.98, np.max(self.density_at_15km) * 1.02)
        ax.set_title('Density over time at h = 15 km (for %s months)' % (self.epoch_text))
        xspan = ax.get_xlim()[1] - ax.get_xlim()[0]
        yspan = ax.get_ylim()[1] - ax.get_ylim()[0]

        for y in np.unique(self.year):
            mjd_start_year = date2mjd(y, 1, 1, 0)
            mjd_half_year = date2mjd(y, 7, 1, 0)
            year_plot = y

            if 1 in np.unique(self.month):
                ax.vlines(mjd_start_year, ax.get_ylim()[0] + 0.12 * yspan, ax.get_ylim()[1], color='0.7', linewidth=1.,
                          linestyles='dotted')
                ax.text(mjd_start_year - xspan * 0.01, ax.get_ylim()[0] + 0.095 * yspan, str(year_plot), rotation=90,
                        color='0.7')

            if 7 in np.unique(self.month):
                ax.vlines(mjd_half_year, ax.get_ylim()[0] + 0.12 * yspan, ax.get_ylim()[1], color='0.7', linewidth=1.,
                          linestyles='dotted')
                ax.text(mjd_half_year - xspan * 0.01, ax.get_ylim()[0] + 0.095 * yspan, str(year_plot + 0.5),
                        rotation=90, color='0.7')

        ax.hlines(np.mean(self.density_at_15km), ax.get_xlim()[0], ax.get_xlim()[1], colors='#336699',
                  linestyle='solid', lw=1., zorder=10)

        fig.savefig(self.output_plot_name + '_at_15_km.eps', bbox_inches='tight')
        fig.savefig(self.output_plot_name + '_at_15_km.png', bbox_inches='tight', dpi=300)

    def plot_differences_wrt_magic(self):
        """
        Function that plots the difference between density models w.r.t the MAGIC Winter model
        :return: 
        """

        print('Computing the averages, std dev and peak-to-peak values for the differences wrt MAGIC Winter model:')
        print('plotting averaged data values for selected epoch')
        fig = plt.figure()
        ax = fig.add_subplot(111)

        eb2 = ax.errorbar(self.x + 175., self.diff_MAGIC[0], yerr=[self.diff_MAGIC[3],
                                                                         self.diff_MAGIC[2]], fmt='o', color='b',
                          capsize=0.5, mec='b', ms=1., label=self.data_server)
        eb2[-1][0].set_linestyle(':')
        ax.errorbar(self.x + 175., self.diff_MAGIC[0], yerr=self.diff_MAGIC[1], fmt=':', color='b',
                    elinewidth=3.1)

        ax.set_title('Relative Difference w.r.t MAGIC W model, for %s months' % (self.epoch_text))
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
        print('Computing the averages, std dev and peak-to-peak values for the differences wrt other model:')
        print('plotting averaged data values for selected epoch')
        if model == 'MW':
            diff = self.diff_MAGIC[0]
            ediff_pp = [self.diff_MAGIC[3], self.diff_MAGIC[2]]
            ediff = self.diff_MAGIC[1]
        elif model == 'PROD3':
            diff = self.diff_PROD3[0]
            ediff_pp = [self.diff_PROD3[3], self.diff_PROD3[2]]
            ediff = self.diff_PROD3[1]
        else:
            print('Wrong model name. It must be "MW" for MAGIC Winter model, or "PROD3" for Paranal model. \n Exiting!')
            sys.exit()

        fig = plt.figure()
        ax = fig.add_subplot(111)

        eb2 = ax.errorbar(self.x + 175., diff, yerr=ediff_pp, fmt='o', color='b', capsize=0.5, mec='b',
                          ms=1., label=self.data_server)
        eb2[-1][0].set_linestyle(':')
        ax.errorbar(self.x + 175., diff, yerr=ediff, fmt=':', color='b',
                    elinewidth=3.1)

        ax.set_title('Relative Difference w.r.t %s model, epoch: %s' % (model,
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

    def plot_models_comparison(self, model=None):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        eb2 = ax.errorbar(self.x, self.averages[0], yerr=[self.averages[3], self.averages[2]],
                          fmt='o', color='b', capsize=0.5, mec='b', ms=1., label=self.data_server)
        eb2[-1][0].set_linestyle(':')
        ax.errorbar(self.x, self.averages[0], yerr=self.averages[1], fmt=':', color='b',
                    elinewidth=3.)

        if model == 'MW':
            ax.plot(self.heightmw * 1000., self.n_mw * np.exp(self.heightmw * 1000. / self.Hs), '-', color='grey',
                    label='MAGIC W')
        elif model == 'PROD3':
            self._get_prod3sim_data()
            ax.plot(self.hprod3 * 1000., self.n_prod3 * np.exp(self.hprod3 * 1000. / self.Hs), '-', color='0.2',
                    label='Prod3 ' + self.observatory)
        elif model == 'both':
            ax.plot(self.heightmw * 1000., self.n_mw * np.exp(self.heightmw * 1000. / self.Hs), '-', color='grey',
                    label='MAGIC W')
            self._get_prod3sim_data()
            ax.plot(self.hprod3 * 1000., self.n_prod3 * np.exp(self.hprod3 * 1000. / self.Hs), '-', color='0.2',
                    label='Prod3 ' + self.observatory)
        else:
            print('Wrong model. It must be "MW" for MAGIC Winter model, or "PROD3" for Paranal model or "both". '
                  '\n Exiting!')
            sys.exit()

        ax.set_title(self.data_server + ' ' + ' ' + str(self.epoch_text))
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
        print('# x[i], avg, rms, p2p_m, p2p_p', file=textfile)
        for i in np.arange(len(self.x)):
            print(self.x[i], self.averages[0][i], self.averages[1][i], self.averages[3][i],
                  self.averages[2][i], self.diff_MAGIC[0][i], self.diff_MAGIC[1][i],
                  self.diff_MAGIC[3][i], self.diff_MAGIC[2][i], file=textfile)
        textfile.close()
