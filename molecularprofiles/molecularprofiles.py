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
import pandas as pd

settings()


class MolecularProfile:
    def __init__(self, data_file, data_server='ECMWF', tag_name='myplots', observatory='north'):

        """
        This class provides with a series of functions to analyze the quality of the data for both
        CTA-North and CTA-South.

        :param data_file: txt file containing the data (string)
        :param tag_name: name to be given to the output files (string)
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

        # P = np.array([1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700, 650, 600, 550, 500, 450, 400,
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

        Input: epoch: (str) can be "winter", "summer", "intermediate", "all". Default is "all"

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
            self.Temp (K)
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

        self.output_plot_name = self.tag_name + '_' + epoch
        self.epoch = epoch

        self.dataframe = pd.read_table(self.data_file.split('.')[0] + '.txt', delimiter=' ')

        if epoch != 'all':
            self.dataframe = select_dataframe_epoch(self.dataframe, epoch)

        self.group_by_p = self.dataframe.groupby('P')
        self.dataframe['n_exp'] = self.dataframe.n /self.Ns * np.exp(self.dataframe.h /self.Hs)
        self.h_avgs = avg_std_dataframe(self.group_by_p, 'h')
        self.n_exp_avgs = avg_std_dataframe(self.group_by_p, 'n_exp')
        self.x = np.linspace(2200., 25000., num=15, endpoint=True)

    def _interpolate_param_to_h(self, param, height):

        interpolated_param = []
        self.group_mjd = self.dataframe.groupby('MJD')

        print("Computing the extrapolation of the values of density:")
        print("(This is to make it easier to compare ECMWF and GDAS, or any other")
        print("weather model)")
        pbar = tqdm(total=len(np.unique(self.dataframe.MJD)))

        for mjd in np.unique(self.dataframe.MJD):
            pbar.update(1)
            h_at_mjd = self.group_mjd.get_group(mjd)['h'].tolist()
            param_at_mjd = self.group_mjd.get_group(mjd)[param].tolist()
            func = interp1d(h_at_mjd, param_at_mjd, kind='cubic', fill_value='extrapolate')

            if type(height) == int or type(height) == float:
                interpolated_param.append(np.float(func(height)))
            else:
                interpolated_param.append(func(height))
        pbar.close()

        print('\n')
        interpolated_param = np.asarray(interpolated_param)
        if type(height) != float:
            interpolated_param_avgs = compute_averages_std(interpolated_param)
            return interpolated_param, interpolated_param_avgs[0], interpolated_param_avgs[1], \
                   interpolated_param_avgs[2], interpolated_param_avgs[3]
        elif type(height) == float or type(height) == int:
            return interpolated_param


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

        pbar = tqdm(total=len(self.dataframe.P))
        for i in np.arange(len(self.dataframe.P)):
            pbar.update(1)
            if air == 'moist':
                Xw.append(MolarFractionWaterVapor(self.dataframe.P[i], self.dataframe.Temp[i], self.dataframe.RH[i]))
                Z.append(Compressibility(self.dataframe.P[i], self.dataframe.Temp[i], Xw[i]))
                rho.append(DensityMoistAir(self.dataframe.P[i] * 100., self.dataframe.Temp[i], Z[i], Xw[i], C))

            elif air == 'dry':
                Z.append(Compressibility(self.dataframe.P[i], self.dataframe.Temp[i], 0.0))
                rho.append(DensityMoistAir(self.dataframe.P[i]*100., self.dataframe.Temp[i], Z[i], 0.0, C))
            else:
                print('Wrong air condition. It must be "moist" or "dry". Aborting!')
                sys.exit()
        pbar.close()

        self.dataframe['n_mass_'+air] = rho
        raw_rho = np.asarray(rho)
        interpolated_rho, rho_avg, e_rho, rho_pp, rho_pm = self._interpolate_param_to_h('n_mass_'+air, self.x)
        return rho, e_rho, rho_pp, rho_pm, raw_rho

    def compute_diff_wrt_model(self, interpolated_density):

        diff_with_magic = []
        diff_with_prod3 = []

        self._get_prod3sim_data()
        x = np.linspace(1000., 25000., num=15, endpoint=True)

        print("Computing the differences of the values of density:")
        for i in tqdm(np.arange(len(interpolated_density))):
            # Percentage counter bar
            # sys.stdout.write('\r')
            # k = 100
            # sys.stdout.write("[%-100s] %d%%" % ('=' * k, k))
            # sys.stdout.flush()
            # ---------------------------
            diff_with_magic.append((interpolated_density[i] - self.func_magic(x))
                                         / interpolated_density[i])
            diff_with_prod3.append((interpolated_density[i] - self.func_prod3(x))
                                         / interpolated_density[i])

        self.diff_with_magic = np.asarray(diff_with_magic)
        self.diff_with_prod3 = np.asarray(diff_with_prod3)
        self.diff_MAGIC = compute_averages_std(self.diff_with_magic)
        self.diff_PROD3 = compute_averages_std(self.diff_with_prod3)

    # =======================================================================================================

    def plot_moist_dry_comparison(self):

        fig, axs = plt.subplots(2,1,sharex=True)
        plt.subplots_adjust(hspace=0)

        av_heights = self.h_avgs
        raw_rhod = self.compute_mass_density(air='dry')[4]
        raw_rhow = self.compute_mass_density(air='moist')[4]

        rel_dif = (raw_rhod - raw_rhow) * 100. / raw_rhod
        self.dataframe['rel_dif_n_mass'] = rel_dif
        self.group_by_p = self.dataframe.groupby('P')
        a_rhod, e_rhod, pp_rhod, pm_rhod = avg_std_dataframe(self.group_by_p, 'n_mass_dry')
        a_rhow, e_rhow, pp_rhow, pm_rhow = avg_std_dataframe(self.group_by_p, 'n_mass_moist')
        a_rel_dif, e_rel_dif, pp_rel_dif, pm_rel_dif = avg_std_dataframe(self.group_by_p,'rel_dif_n_mass')

        axs[0].errorbar(av_heights[0], a_rhod, yerr=e_rhod, fmt=':', color='#ff7f0e', elinewidth=3)
        axs[0].errorbar(av_heights[0], a_rhow, yerr=e_rhow, fmt=':', color='#1f77b4', elinewidth=3)
        ebw = axs[0].errorbar(av_heights[0], a_rhow, yerr=[pm_rhow, pp_rhow], fmt='o', color='#1f77b4', capsize=0.5, mec='#1f77b4', ms=2., label='$\\rho_w$ (moist air)')
        ebd = axs[0].errorbar(av_heights[0], a_rhod, yerr=[pm_rhod, pp_rhod], fmt='o', color='#ff7f0e', capsize=0.5, mec='#ff7f0e', ms=2., label='$\\rho_d$ (dry air)')
        ebd[-1][0].set_linestyle(':')
        ebw[-1][0].set_linestyle(':')

        axs[0].legend(loc='best')
        #axs[0].set_ylabel('$\\rho$ * exp(h/H$_{\\rm s}$) [kg m$^{-3}$]')
        axs[0].axes.tick_params(direction='in')

        axs[1].axes.tick_params(direction='inout', top='on')
        axs[1].set_xlabel('height [m]')
        axs[1].set_yscale('log')
        axs[1].errorbar(av_heights[0], a_rel_dif, yerr=e_rel_dif, color='#2ca02c', ms=2.)
        ebrd = axs[1].errorbar(av_heights[0], a_rel_dif, yerr=[pm_rel_dif, pp_rel_dif], capsize=0.5, ms=2.)
        ebrd[-1][0].set_linestyle(':')
        axs[1].set_ylabel('rel. diff [%]')
        axs[1].set_ylim(1.e-4, 5.e-1)

        fig.savefig('dry_vs_moist_air_density_RH_lt_0.80_h_gt_2200_rel_dif_error.png', bbox_inches='tight', dpi=300)
        plt.show()

    def plot_moist_dry_comparison_interpolated(self):

        fig, axs = plt.subplots(2, 1, sharex=True)
        plt.subplots_adjust(hspace=0)

        av_heights = self.x
        a_rhod, e_rhod, pp_rhod, pm_rhod, raw_rhod = self.compute_mass_density(air='dry')
        a_rhow, e_rhow, pp_rhow, pm_rhow, raw_rhow = self.compute_mass_density(air='moist')

        rel_dif = (a_rhod - a_rhow) * 100. / a_rhod

        axs[0].errorbar(av_heights, a_rhod, yerr=e_rhod, fmt=':', color='#ff7f0e', elinewidth=3)
        axs[0].errorbar(av_heights, a_rhow, yerr=e_rhow, fmt=':', color='#1f77b4', elinewidth=3)
        ebw = axs[0].errorbar(av_heights, a_rhow, yerr=[pm_rhow, pp_rhow], fmt='o', color='#1f77b4', capsize=0.5,
                              mec='#1f77b4', ms=2., label='$\\rho_w$ (moist air)')
        ebd = axs[0].errorbar(av_heights, a_rhod, yerr=[pm_rhod, pp_rhod], fmt='o', color='#ff7f0e', capsize=0.5,
                              mec='#ff7f0e', ms=2., label='$\\rho_d$ (dry air)')
        ebd[-1][0].set_linestyle(':')
        ebw[-1][0].set_linestyle(':')

        axs[0].legend(loc='best')
        axs[0].set_ylabel('$\\rho$ * exp(h/H$_{\\rm s}$) [kg m$^{-3}$]')
        axs[0].axes.tick_params(direction='in')

        axs[1].axes.tick_params(direction='inout', top='on')
        axs[1].set_xlabel('height [m]')
        axs[1].set_yscale('log')
        axs[1].errorbar(av_heights, rel_dif, color='#2ca02c', ms=2.)
        ebrd = axs[1].errorbar(av_heights, rel_dif, capsize=0.5, ms=2.)
        ebrd[-1][0].set_linestyle(':')
        axs[1].set_ylabel('rel. diff [\%]')
        axs[1].set_ylim(1.e-4, 5.e-1)

        fig.savefig('interpolated_dry_vs_moist_air_density_RH_lt_0.80_h_gt_2200_rel_dif_error.png', bbox_inches='tight', dpi=300)
        plt.show()

    def plot_average_at_15km(self):
        """
        Function that produces a plot of the averaged density at 15 km
        :return: 
        """

        print('plotting data at 15 km in 1-day bins:')

        fig = plt.figure()
        ax = fig.add_subplot(111)
        density_at_15km = self._interpolate_param_to_h('n_exp', 15000.)
        ax.plot(np.unique(self.dataframe.MJD), density_at_15km, 'o', color='#99CCFF', markersize=1.2,
                label=self.data_server + ' ' + self.observatory, alpha=0.8)

        ax.legend(loc='best', numpoints=1)
        ax.set_xlabel('MJD')
        ax.set_ylabel('$n_{\\rm day}/N_{\\rm s} \\cdot e^{(h/H_{\\rm s})}$')
        ax.set_ylim(np.min(density_at_15km) * 0.98, np.max(density_at_15km) * 1.02)
        ax.set_title('Density over time at h = 15 km (for %s months)' % (self.epoch))
        xspan = ax.get_xlim()[1] - ax.get_xlim()[0]
        yspan = ax.get_ylim()[1] - ax.get_ylim()[0]

        for y in np.unique(self.dataframe.year):
            mjd_start_year = date2mjd(y, 1, 1, 0)
            mjd_half_year = date2mjd(y, 7, 1, 0)
            year_plot = y

            if 1 in np.unique(self.dataframe.month):
                ax.vlines(mjd_start_year, ax.get_ylim()[0] + 0.12 * yspan, ax.get_ylim()[1], color='0.7', linewidth=1.,
                          linestyles='dotted')
                ax.text(mjd_start_year - xspan * 0.01, ax.get_ylim()[0] + 0.095 * yspan, str(year_plot), rotation=90,
                        color='0.7')

            if 7 in np.unique(self.dataframe.month):
                ax.vlines(mjd_half_year, ax.get_ylim()[0] + 0.12 * yspan, ax.get_ylim()[1], color='0.7', linewidth=1.,
                          linestyles='dotted')
                ax.text(mjd_half_year - xspan * 0.01, ax.get_ylim()[0] + 0.095 * yspan, str(year_plot + 0.5),
                        rotation=90, color='0.7')

        ax.hlines(np.mean(density_at_15km), ax.get_xlim()[0], ax.get_xlim()[1], colors='#336699',
                  linestyle='solid', lw=1., zorder=10)

        fig.savefig(self.output_plot_name + '_at_15_km.eps', bbox_inches='tight')
        fig.savefig(self.output_plot_name + '_at_15_km.png', bbox_inches='tight', dpi=300)

    def plot_differences_wrt_magic(self):
        """
        Function that plots the difference between density models w.r.t the MAGIC Winter model
        :return: 
        """
        print('plotting averaged data values for selected epoch')
        fig = plt.figure()
        ax = fig.add_subplot(111)

        eb2 = ax.errorbar(self.x + 175., self.diff_MAGIC[0], yerr=[self.diff_MAGIC[3],
                                                                         self.diff_MAGIC[2]], fmt='o', color='b',
                          capsize=0.5, mec='b', ms=1., label=self.data_server)
        eb2[-1][0].set_linestyle(':')
        ax.errorbar(self.x + 175., self.diff_MAGIC[0], yerr=self.diff_MAGIC[1], fmt=':', color='b',
                    elinewidth=3.1)

        ax.set_title('Relative Difference w.r.t MAGIC W model, for %s months' % (self.epoch))
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

    def plot_differences_wrt_model(self, model):
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

        ax.set_title('Relative Difference w.r.t %s model, epoch: %s' % (model, self.epoch))
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

    def plot_models_comparison(self, model=None, interpolate=False):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if interpolate:
            raw_n_exp, avg_n_exp, e_n_exp, pp_n_exp, pm_n_exp = self._interpolate_param_to_h('n_exp', self.x)
            eb2 = ax.errorbar(self.x, avg_n_exp, yerr=[pm_n_exp, pp_n_exp],
                          fmt='o', color='b', capsize=0.5, mec='b', ms=1., label=self.data_server)
            eb2[-1][0].set_linestyle(':')
            ax.errorbar(self.x, avg_n_exp, yerr=e_n_exp, fmt=':', color='b', elinewidth=3., label=None)
        else:
            eb2 = ax.errorbar(self.h_avgs[0], self.n_exp_avgs[0], xerr=self.h_avgs[1], yerr=[self.n_exp_avgs[3],
                                                                                             self.n_exp_avgs[2]],
                              fmt='o', color='b', capsize=0.5, mec='b', ms=1., label=self.data_server)
            eb2[-1][0].set_linestyle(':')
            ax.errorbar(self.h_avgs[0], self.n_exp_avgs[0], xerr=self.h_avgs[1], yerr=self.n_exp_avgs[1], fmt=':',
                        color='b', elinewidth=3., label=None)

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

        ax.set_title(self.data_server + ' ' + ' ' + str(self.epoch))
        ax.set_xlabel('h a.s.l. [m]')
        ax.set_ylabel('$n_{\\rm day}/N_{\\rm s} \\cdot e^{(h/H_{\\rm s})}$')
        ax.set_xlim(0., 25100.)
        ax.set_ylim(0.4, 1.2)
        ax.xaxis.set_minor_locator(MultipleLocator(1000))
        ax.xaxis.set_major_locator(MultipleLocator(2000))
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
        ax.legend(loc='best', numpoints=1)
        ax.grid(which='both', axis='y', color='0.8')
#        fig.savefig('comparison_' + self.output_plot_name + '.eps', bbox_inches='tight')
        fig.savefig('comparison_' + self.output_plot_name + '.png', bbox_inches='tight', dpi=300)

    def print_to_text_file(self):
        textfile = open(self.output_plot_name + '_to_text_file.txt', 'w')
        print('# x[i], avg, rms, p2p_m, p2p_p', file=textfile)
        for i in np.arange(len(self.x)):
            print(self.x[i], self.averages[0][i], self.averages[1][i], self.averages[3][i],
                  self.averages[2][i], self.diff_MAGIC[0][i], self.diff_MAGIC[1][i],
                  self.diff_MAGIC[3][i], self.diff_MAGIC[2][i], file=textfile)
        textfile.close()
