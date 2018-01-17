import os
import os.path
import sys
import matplotlib.pyplot as plt
from tqdm import *
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import interp1d
from molecularprofiles.utils.grib_utils import *
from molecularprofiles.utils.read_txtfile import *
from molecularprofiles.plot_settings import settings
from molecularprofiles.utils.magic_winter_profile import heightmw, rhomw
from molecularprofiles.utils.meteorological_constants import *

settings()


class EcmwfMolecularProfile:
    def __init__(self, file_ecmwf, tag_name='myplots', epoch_text='all',
                 label_ecmwf='ECMWF', observatory='north'):

        """
        This class provides with a series of functions to analyze the quality of the ECMWF data for both
        CTA-North and CTA-South.
        
        :param file_ecmwf: txt file containing the ECMWF data (string) 
        :param tag_name: name to be given to the output files (string)
        :param epoch_text: season to be analyzed. Valid values are: "all", "winter", "summer", "intermediate" (string) 
        :param label_ecmwf: label to be put in some of the output plots (string)
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

        self.file_ecmwf = file_ecmwf
        self.epoch_text = epoch_text
        self.tag_name = tag_name
        self.label_ecmwf = label_ecmwf
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

        # TODO: change the directory definition. Make an import since everything is now recognised as package

        MOLECULARPROFILES_DIR = os.environ.get('MOLECULARPROFILES_DIR')
        ECMWF_DIR = MOLECULARPROFILES_DIR + 'ecmwf_scipts/'
        PROD3_DIR = MOLECULARPROFILES_DIR + 'Prod3b_simulations/'
        # Prod3 Simulations (based on NRLMSISE)
        if self.observatory == 'north':
            prod3 = open(PROD3_DIR+'/atmprof36_lapalma.dat')
            self.hprod3, nprod3, thickprod3, n1prod3, tempprod3, pprod3 = np.loadtxt(prod3, usecols=(0,1,2,3,4,5),
                                                                                     unpack=True)
        elif self.observatory == 'shouth':
            prod3 = open(PROD3_DIR+'/atmprof26_paranal.dat')
            self.hprod3, nprod3, thickprod3, n1prod3, tempprod3, pprod3 = np.loadtxt(prod3, usecols=(0, 1, 2, 3, 4, 5),
                                                                                     unpack=True)
        else:
            print ('WRONG observatory. It must be "north" or "south" ')
            raise SystemExit

        self.n_prod3 = nprod3 * 816.347
        self.func_prod3 = interp1d(self.hprod3 * 1000., self.n_prod3 * np.exp(self.hprod3 * 1000. / self.Hs),
                                   kind='cubic')


    def get_ecmwf_data(self):

        """
        Function that reads ECMWF txt input files and returns quantities ready to plot.
        If the input ecmwf filename does not exist, the program searches for the same
        input name with .grib extension and extracts the data from it
        :return: ()
        """

        if not os.path.exists((self.file_ecmwf).split('.')[0]+'.txt'):
            grib_file = self.file_ecmwf
            readgribfile2text(grib_file, self.observatory, gridstep=0.75)

        self.mjd_ecmwf, self.year_ecmwf, self.month_ecmwf, self.day_ecmwf, self.hour_ecmwf, self.p_ecmwf, self.h_ecmwf, \
        self.n_ecmwf = read_file(self.file_ecmwf.split('.')[0]+'.txt', self.epoch_text)

        interpolated_density_ecmwf = []
        ecmwf_density_at_15km = []
        mjd_at_15km = []
        month_at_15km = []
        mjd = self.mjd_ecmwf[0]
        self.x = np.linspace(1000., 25000., num=15, endpoint=True)
        step_hours = self.mjd_ecmwf[37] - self.mjd_ecmwf[0]
        steps = (np.max(self.mjd_ecmwf) - mjd) / step_hours

        pbar = tqdm(total = steps + 1)

        print("Computing the extrapolation of the values of density for ECMWF:")
        while mjd < (np.max(self.mjd_ecmwf) + 0.25):
            # Percentage counter bar
            pbar.update(1)
            # ---------------------------
            func_ecmwf = interp1d(self.h_ecmwf[self.mjd_ecmwf == mjd], self.n_ecmwf[self.mjd_ecmwf == mjd] / self.Ns *
                                  np.exp(self.h_ecmwf[self.mjd_ecmwf == mjd] / self.Hs), kind='cubic')

            int_dens_ecwmf = func_ecmwf(self.x)
            interpolated_density_ecmwf.append(int_dens_ecwmf)
            ecmwf_density_at_15km.append(func_ecmwf(self.x[8]))
            mjd_at_15km.append(mjd)
            month_at_15km.append(self.month_ecmwf[self.mjd_ecmwf == mjd])

            mjd += step_hours
        pbar.close()

        print('\n')
        self.interpolated_density_ecmwf = np.asarray(interpolated_density_ecmwf)
        self.ecmwf_density_at_15km = np.asarray(ecmwf_density_at_15km)
        self.mjd_at_15km_ecmwf = np.asarray(mjd_at_15km)
        self.ecmwf_averages = compute_averages_std(self.interpolated_density_ecmwf)

    def compute_diff_wrt_model(self):

        diff_ecmwf_with_magic = []
        diff_ecmwf_with_prod3 = []

        self._get_prod3sim_data()
        self.x = np.linspace(1000., 25000., num=15, endpoint=True)

        print("Computing the differences of the values of density for ECMWF:")
        for i in tqdm(np.arange(len(self.interpolated_density_ecmwf))):
            # Percentage counter bar
            #sys.stdout.write('\r')
            #k = 100
            #sys.stdout.write("[%-100s] %d%%" % ('=' * k, k))
            #sys.stdout.flush()
            # ---------------------------
            diff_ecmwf_with_magic.append((self.interpolated_density_ecmwf[i] - self.func_magic(self.x))
                                         / self.interpolated_density_ecmwf[i])
            diff_ecmwf_with_prod3.append((self.interpolated_density_ecmwf[i] - self.func_prod3(self.x))
                                         / self.interpolated_density_ecmwf[i])

        self.diff_ecmwf_with_magic = np.asarray(diff_ecmwf_with_magic)
        self.diff_ecmwf_with_prod3 = np.asarray(diff_ecmwf_with_prod3)
        self.diff_ecmwf_MAGIC = compute_averages_std(self.diff_ecmwf_with_magic)
        self.diff_ecmwf_PROD3 = compute_averages_std(self.diff_ecmwf_with_prod3)

    # =======================================================================================================

    def plot_average_at_15km(self):
        """
        Function that produces a plot of the averaged density at 15 km for ECMWF data
        :return: 
        """

        print('plotting data at 15 km in 1-day bins:')

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.mjd_at_15km_ecmwf, self.ecmwf_density_at_15km, 'o', color='#99CCFF', markersize=1.2,
                label=self.label_ecmwf+' '+self.observatory, alpha=0.8)

        ax.legend(loc='best', numpoints=1)
        ax.set_xlabel('MJD')
        ax.set_ylabel('$n_{\\rm day}/N_{\\rm s} \\cdot e^{(h/H_{\\rm s})}$')
        ax.set_ylim(np.min(self.ecmwf_density_at_15km)*0.98, np.max(self.ecmwf_density_at_15km)*1.02)
        ax.set_xlim(np.min(self.mjd_at_15km_ecmwf)- 10, np.max(self.mjd_at_15km_ecmwf) + 10)
        ax.set_title('Density over time at h = 15 km (%s)' %(self.epoch_text))
        xspan = ax.get_xlim()[1] - ax.get_xlim()[0]

        for y in np.unique(self.year_ecmwf):
            mjd_start_year = date2mjd(y, 1, 1, 0)
            mjd_half_year = date2mjd(y, 7, 1, 0)
            year_plot = y

            if 1 in np.unique(self.month_ecmwf):
                ax.vlines(mjd_start_year, ax.get_ylim()[0]*1.02, ax.get_ylim()[1], color='0.7', linewidth=1.,
                          linestyles='dotted')
                ax.text(mjd_start_year - xspan*0.01, ax.get_ylim()[0]*1.015, str(year_plot), rotation=90, color='0.7')

            if 7 in np.unique(self.month_ecmwf):
                ax.vlines(mjd_half_year,  ax.get_ylim()[0], ax.get_ylim()[1], color='0.7', linewidth=1.,
                          linestyles='dotted')
                ax.text(mjd_half_year, 0.77, str(year_plot), rotation=90, color='0.7')

        ax.hlines(np.mean(self.ecmwf_density_at_15km), ax.get_xlim()[0], ax.get_xlim()[1], colors='#336699',
                  linestyle='solid', lw=1., zorder=10)

        fig.savefig(self.output_plot_name + '_at_15_km.eps', bbox_inches='tight')
        fig.savefig(self.output_plot_name + '_at_15_km.png', bbox_inches='tight', dpi=300)

    def plot_differences_wrt_magic(self):
        """
        Function that plots the difference between the ECMWF density models w.r.t the MAGIC Winter model
        :return: 
        """

        print('Computing the averages, std dev and peak-to-peak values for the differences wrt MAGIC Winter model:')
        print('plotting averaged data values for selected epoch')
        fig = plt.figure()
        ax = fig.add_subplot(111)

        eb2 = ax.errorbar(self.x + 175., self.diff_ecmwf_MAGIC[0], yerr=[self.diff_ecmwf_MAGIC[3],
                                                                         self.diff_ecmwf_MAGIC[2]], fmt='o', color='b',
                          capsize=0.5, mec='b', ms=1., label=self.label_ecmwf)
        eb2[-1][0].set_linestyle(':')
        ax.errorbar(self.x + 175., self.diff_ecmwf_MAGIC[0], yerr=self.diff_ecmwf_MAGIC[1], fmt=':', color='b',
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
            diff_ecmwf = self.diff_ecmwf_MAGIC[0]
            ediff_ecmwf_pp = [self.diff_ecmwf_MAGIC[3], self.diff_ecmwf_MAGIC[2]]
            ediff_ecmwf = self.diff_ecmwf_MAGIC[1]
        elif model == 'PROD3':
            diff_ecmwf = self.diff_ecmwf_PROD3[0]
            ediff_ecmwf_pp = [self.diff_ecmwf_PROD3[3], self.diff_ecmwf_PROD3[2]]
            ediff_ecmwf = self.diff_ecmwf_PROD3[1]
        else:
            print('Wrong model name. It must be "MW" for MAGIC Winter model, or "PROD3" for Paranal model. \n Exiting!')
            sys.exit()

        fig = plt.figure()
        ax = fig.add_subplot(111)

        eb2 = ax.errorbar(self.x + 175., diff_ecmwf, yerr=ediff_ecmwf_pp, fmt='o', color='b', capsize=0.5, mec='b',
                          ms=1., label=self.label_ecmwf)
        eb2[-1][0].set_linestyle(':')
        ax.errorbar(self.x + 175., diff_ecmwf, yerr=ediff_ecmwf, fmt=':', color='b',
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
        else:
            print('Wrong model. It must be "MW" for MAGIC Winter model, or "PROD3" for Paranal model or "both". '
                  '\n Exiting!')
            sys.exit()

        ax.set_title(self.label_ecmwf + ' ' + ' ' + str(self.epoch_text))
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
        print('# x[i], ecmwf_avg, ecmwf_rms, ecmwf_p2p_m, ecmwf_p2p_p', file=textfile)
        for i in np.arange(len(self.x)):
            print(self.x[i], self.ecmwf_averages[0][i], self.ecmwf_averages[1][i], self.ecmwf_averages[3][i],
                  self.ecmwf_averages[2][i], self.diff_ecmwf_MAGIC[0][i], self.diff_ecmwf_MAGIC[1][i],
                  self.diff_ecmwf_MAGIC[3][i], self.diff_ecmwf_MAGIC[2][i], file=textfile)
        textfile.close()
