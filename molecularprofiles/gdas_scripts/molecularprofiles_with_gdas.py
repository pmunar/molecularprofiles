import os
import os.path
import sys
import matplotlib.pyplot as plt
import time
from tqdm import *
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import interp1d
from molecularprofiles.utils.grib_utils import *
from molecularprofiles.utils.read_txt_gribfile import *
from molecularprofiles.plot_settings import settings
from molecularprofiles.utils.magic_winter_profile import heightmw, rhomw
from molecularprofiles.utils.meteorological_constants import *

settings()


class gdasMolecularProfile:
    def __init__(self, file_gdas, tag_name='myplots', epoch_text='all',
                 label_gdas='gdas', observatory='north'):

        """
        This class provides with a series of functions to analyze the quality of the gdas data for both
        CTA-North and CTA-South.

        :param file_gdas: txt file containing the gdas data (string) 
        :param tag_name: name to be given to the output files (string)
        :param epoch_text: season to be analyzed. Valid values are: "all", "winter", "summer", "intermediate" (string) 
        :param label_gdas: label to be put in some of the output plots (string)
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

        self.file_gdas = file_gdas
        self.epoch_text = epoch_text
        self.tag_name = tag_name
        self.label_gdas = label_gdas
        self.observatory = observatory
        self.output_plot_name = self.tag_name + '_' + self.epoch_text

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

        MOLECULARPROFILES_DIR = os.environ.get('MOLECULARPROFILES_DIR')
        gdas_DIR = MOLECULARPROFILES_DIR + '/molecularprofiles/gdas_scipts/'
        PROD3_DIR = MOLECULARPROFILES_DIR + '/molecularprofiles/Prod3b_simulations/'
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

    def get_gdas_data(self):

        """
        Function that reads gdas txt input files and returns quantities ready to plot.
        If the input gdas filename does not exist, the program searches for the same
        input name with .grib extension and extracts the data from it
        :return: ()
        """

        if not os.path.exists((self.file_gdas).split('.')[0] + '.txt'):
            grib_file = self.file_gdas
            readgribfile2text(grib_file, self.observatory, gridstep=1.0)

        self.mjd_gdas, self.year_gdas, self.month_gdas, self.day_gdas, self.hour_gdas, self.p_gdas, self.h_gdas, \
        self.n_gdas = read_file(self.file_gdas.split('.')[0] + '.txt', self.epoch_text)

        interpolated_density_gdas = []
        gdas_density_at_15km = []
        mjd_at_15km = []
        self.x = np.linspace(1000., 25000., num=15, endpoint=True)

        print("Computing the extrapolation of the values of density for gdas:")
        pbar = tqdm(total = len(np.unique(self.mjd_gdas)))

        for mjd in np.unique(self.mjd_gdas):
            pbar.update(1)
            func_gdas = interp1d(self.h_gdas[self.mjd_gdas == mjd], self.n_gdas[self.mjd_gdas == mjd] / self.Ns
                                   * np.exp(self.h_gdas[self.mjd_gdas == mjd] / self.Hs), kind='cubic')

            int_dens_ecwmf = func_gdas(self.x)
            interpolated_density_gdas.append(int_dens_ecwmf)
            gdas_density_at_15km.append(func_gdas(self.x[8]))
            mjd_at_15km.append(mjd)
        pbar.close()

        print('\n')
        self.interpolated_density_gdas = np.asarray(interpolated_density_gdas)
        self.gdas_density_at_15km = np.asarray(gdas_density_at_15km)
        self.mjd_at_15km_gdas = np.asarray(mjd_at_15km)
        self.gdas_averages = compute_averages_std(self.interpolated_density_gdas)

    def compute_diff_wrt_model(self):

        diff_gdas_with_magic = []
        diff_gdas_with_prod3 = []

        self._get_prod3sim_data()
        self.x = np.linspace(1000., 25000., num=15, endpoint=True)

        print("Computing the differences of the values of density for gdas:")
        for i in tqdm(np.arange(len(self.interpolated_density_gdas))):
            # Percentage counter bar
            # sys.stdout.write('\r')
            # k = 100
            # sys.stdout.write("[%-100s] %d%%" % ('=' * k, k))
            # sys.stdout.flush()
            # ---------------------------
            diff_gdas_with_magic.append((self.interpolated_density_gdas[i] - self.func_magic(self.x))
                                         / self.interpolated_density_gdas[i])
            diff_gdas_with_prod3.append((self.interpolated_density_gdas[i] - self.func_prod3(self.x))
                                         / self.interpolated_density_gdas[i])

        self.diff_gdas_with_magic = np.asarray(diff_gdas_with_magic)
        self.diff_gdas_with_prod3 = np.asarray(diff_gdas_with_prod3)
        self.diff_gdas_MAGIC = compute_averages_std(self.diff_gdas_with_magic)
        self.diff_gdas_PROD3 = compute_averages_std(self.diff_gdas_with_prod3)

    # =======================================================================================================

    def plot_average_at_15km(self):
        """
        Function that produces a plot of the averaged density at 15 km for gdas data
        :return: 
        """

        print('plotting data at 15 km in 1-day bins:')

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.mjd_at_15km_gdas, self.gdas_density_at_15km, 'o', color='#99CCFF', markersize=0.75,
                label=self.label_gdas + ' ' + self.observatory, alpha=0.7)

        mjd_start_year = np.array([date2mjd(self.year_gdas[0], 1, 1, 0)])
        mjd_half_year = np.array([date2mjd(self.year_gdas[0], 7, 1, 0)])
        year_plot = np.array([self.year_gdas[0]])

        if mjd_start_year.all():
            for i in np.arange(len(mjd_start_year)):
                ax.vlines(mjd_start_year[i], 0.775, 1.0, color='0.7', linewidth=1.)
                ax.text(mjd_start_year[i] - 7.25, 0.79, str(year_plot[i]), rotation=90, color='0.7')
        if mjd_half_year.all():
            for i in np.arange(len(mjd_half_year)):
                ax.vlines(mjd_half_year[i], 0.775, 1.0, color='0.7', linewidth=1.)
                ax.text(mjd_half_year[i] - 7.25, 0.79, str(year_plot[i]), rotation=90, color='0.7')

        ax.hlines(np.mean(self.gdas_density_at_15km), np.min(self.mjd_at_15km_gdas) - 10.,
                  np.max(self.mjd_at_15km_gdas) + 10.,
                  colors='b', linestyle='solid', lw=2., zorder=10)

        ax.legend(loc='best', numpoints=1)
        ax.set_xlabel('MJD')
        ax.set_ylabel('$n_{\\rm day}/N_{\\rm s} \\cdot e^{(h/H_{\\rm s})}$')
        ax.set_ylim(0.757, 0.939)
        ax.set_title('Density over time at h = 15 km (' + str(self.epoch_text) + ')')
        fig.savefig(self.output_plot_name + '_at_15_km.eps', bbox_inches='tight')
        fig.savefig(self.output_plot_name + '_at_15_km.png', bbox_inches='tight', dpi=300)

    def plot_differences_wrt_magic(self):
        """
        Function that plots the difference between the gdas density models w.r.t the MAGIC Winter model
        :return: 
        """

        print('Computing the averages, std dev and peak-to-peak values for the differences wrt MAGIC Winter model:')
        print('plotting averaged data values for selected epoch')
        fig = plt.figure()
        ax = fig.add_subplot(111)

        eb2 = ax.errorbar(self.x + 175., self.diff_gdas_MAGIC[0], yerr=[self.diff_gdas_MAGIC[3],
                                                                         self.diff_gdas_MAGIC[2]], fmt='o', color='b',
                          capsize=0.5, mec='b', ms=1., label=self.label_gdas)
        eb2[-1][0].set_linestyle(':')
        ax.errorbar(self.x + 175., self.diff_gdas_MAGIC[0], yerr=self.diff_gdas_MAGIC[1], fmt=':', color='b',
                    elinewidth=3.1)

        ax.set_title('Relative Difference w.r.t MAGIC W model %s' % (self.epoch_text))
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

    def plot_differences_wrt_model(self, model=None):
        print('Computing the averages, std dev and peak-to-peak values for the differences wrt model:')
        print('plotting averaged data values for selected epoch')
        print('NEED TO FINisH THIS!')
        if model == 'MW':
            diff_gdas = self.diff_gdas_MAGIC[0]
            ediff_gdas_pp = [self.diff_gdas_MAGIC[3], self.diff_gdas_MAGIC[2]]
            ediff_gdas = self.diff_gdas_MAGIC[1]
        elif model == 'PROD3':
            diff_gdas = self.diff_gdas_PROD3[0]
            ediff_gdas_pp = [self.diff_gdas_PROD3[3], self.diff_gdas_PROD3[2]]
            ediff_gdas = self.diff_gdas_PROD3[1]
        else:
            print('Wrong model. It must be "MW" for MAGIC Winter model, or "PROD3" for Paranal model. \n Exiting!')
            sys.exit()

        fig = plt.figure()
        ax = fig.add_subplot(111)

        eb2 = ax.errorbar(self.x + 175., diff_gdas, yerr=ediff_gdas_pp, fmt='o', color='b', capsize=0.5, mec='b',
                          ms=1., label=self.label_gdas)
        eb2[-1][0].set_linestyle(':')
        ax.errorbar(self.x + 175., diff_gdas, yerr=ediff_gdas, fmt=':', color='b',
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
        fig.savefig('differences_wrt_'+str(model)+'_' + self.output_plot_name + '.eps', bbox_inches='tight')
        fig.savefig('differences_wrt_'+str(model)+'_' + self.output_plot_name + '.png', bbox_inches='tight', dpi=300)

    def plot_models_comparison(self, model=None):
        fig = plt.figure()
        ax = fig.add_subplot(111)

        eb2 = ax.errorbar(self.x + 175., self.gdas_averages[0], yerr=[self.gdas_averages[3], self.gdas_averages[2]],
                          fmt='o', color='b', capsize=0.5, mec='b', ms=1., label=self.label_gdas)
        eb2[-1][0].set_linestyle(':')
        ax.errorbar(self.x + 175., self.gdas_averages[0], yerr=self.gdas_averages[1], fmt=':', color='b',
                    elinewidth=3.)

        if model == 'MW':
            ax.plot(self.heightmw * 1000., self.n_mw * np.exp(self.heightmw * 1000. / self.Hs), '-', color='grey',
                    label='MAGIC W')
        elif model == 'PROD3':
            ax.plot(self.hprod3 * 1000., self.n_prod3 * np.exp(self.hprod3 * 1000. / self.Hs), '-', color='0.2',
                    label='Prod3 ' + self.observatory)
        elif model == 'both':
            ax.plot(self.heightmw * 1000., self.n_mw * np.exp(self.heightmw * 1000. / self.Hs), ':', color='0.8',
                    label='MAGIC W')
            ax.plot(self.hprod3 * 1000., self.n_prod3 * np.exp(self.hprod3 * 1000. / self.Hs), '-', color='0.2',
                    label='Prod3 ' + self.observatory)
        else:
            print('Wrong model. It must be "MW" for MAGIC Winter model, or "PROD3" for Paranal model. \n Exiting')
            sys.exit()

        ax.set_title(self.label_gdas + ' ' + ' ' + str(self.epoch_text))
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
        print('# x[i], gdas_avg, gdas_rms, gdas_p2p_m, gdas_p2p_p', file=textfile)
        for i in np.arange(len(self.x)):
            print(self.x[i], self.gdas_averages[0][i], self.gdas_averages[1][i], self.gdas_averages[3][i],
                  self.gdas_averages[2][i], self.diff_gdas_MAGIC[0][i], self.diff_gdas_MAGIC[1][i],
                  self.diff_gdas_MAGIC[3][i], self.diff_gdas_MAGIC[2][i], file=textfile)
        textfile.close()
