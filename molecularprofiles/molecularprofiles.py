#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# This is a class that offers a set of functions to work with meteorological data
# in txt table format, separated with spaces.
#
# It makes use of the pandas dataframe objects in order to be more efficient
#
#
# Created by Pere Munar-Adrover
# email: pere.munaradrover@gmail.com
#
#
# Version 2.0.1
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

import os
import sys
import matplotlib.pyplot as plt
from tqdm import *
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import interp1d
#from molecularprofiles.utils.observatory import *
from molecularprofiles.utils.grib_utils import *
from molecularprofiles.utils.plot_settings import settings
from molecularprofiles.aux.magic_winter_profile import heightmw, rhomw
from molecularprofiles.aux.meteorological_constants import *
from LIDAR_Analysis.humidity import *
from LIDAR_Analysis.rayleigh import Rayleigh
import pandas as pd

settings()


class MolecularProfile:
    def __init__(self, data_file, data_server='server', tag_name='myplots', observatory='north'):

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
        self.data_server = data_server
        self.tag_name = tag_name + '_' + self.data_server
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
        PROD3_DIR = MOLECULARPROFILES_DIR + '/molecularprofiles/aux/Prod3b_simulations'
        # Prod3 Simulations (based on NRLMSISE)
        if self.observatory == 'north':
            prod3 = open(PROD3_DIR + '/atmprof36_lapalma.dat')
            self.hprod3, nprod3, thickprod3, n1prod3, tempprod3, pprod3 = np.loadtxt(prod3, usecols=(0, 1, 2, 3, 4, 5),
                                                                                     unpack=True)
        elif self.observatory == 'south':
            prod3 = open(PROD3_DIR + '/atmprof26_paranal.dat')
            self.hprod3, nprod3, thickprod3, n1prod3, tempprod3, pprod3 = np.loadtxt(prod3, usecols=(0, 1, 2, 3, 4, 5),
                                                                                     unpack=True)
        else:
            print('WRONG observatory. It must be "north" or "south" ')
            raise SystemExit

        self.n_prod3 = nprod3 * 816.347
        self.func_prod3 = interp1d(self.hprod3 * 1000., self.n_prod3 * np.exp(self.hprod3 * 1000. / self.Hs),
                                   kind='cubic')

    def _interpolate_simple(self, x_param, y_param, new_x_param):
        func = interp1d(x_param, y_param, kind='cubic', fill_value='extrapolate')
        return func(new_x_param)

    def _interpolate_param_to_h(self, param, height):

        interpolated_param = []
        group_mjd = self.dataframe.groupby('MJD')

        print("Computing the extrapolation of the values of density:")
        print("(This is to make it easier to compare ECMWF and GDAS, or any other")
        print("weather model)")
        pbar = tqdm(total=len(np.unique(self.dataframe.MJD)))

        for mjd in np.unique(self.dataframe.MJD):
            pbar.update(1)
            h_at_mjd = group_mjd.get_group(mjd)['h'].tolist()
            param_at_mjd = group_mjd.get_group(mjd)[param].tolist()
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

    def get_data(self, epoch='all', years=None, months=None, hours=None, altitude=[], select_good_weather=False,
                 RH_lim=100., W_lim = 10000., epoch_by_density=False, n_exp_minvalue=0., n_exp_maxvalue=1.,
                 epoch_by_density_name='winter'):

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
            if self.data_server == 'GDAS':
                gridstep = 1.0
            elif self.data_server == 'ECMWF':
                gridstep = 0.75
            readgribfile2text(grib_file, gridstep, self.observatory)

        self.output_plot_name = self.tag_name + '_' + epoch
        self.epoch = epoch

        self.dataframe = pd.read_csv(self.data_file.split('.')[0] + '.txt', sep=' ', comment='#')
        self.dataframe['n_exp'] = self.dataframe.n / self.Ns * np.exp(self.dataframe.h / self.Hs)
        # Altitude filtering:
        if altitude != [] and len(altitude) == 2 and altitude[0] < altitude[1]:
            altitude_cond = (self.dataframe.h >= altitude[0]) & (self.dataframe.h < altitude[1])
            self.dataframe = self.dataframe[altitude_cond]
            self.x = np.linspace(altitude[0], altitude[1], num=15, endpoint=True)
        elif altitude != [] and len(altitude) == 1:
            altitude_cond = (self.dataframe.h < altitude[0])
            self.dataframe = self.dataframe[altitude_cond]
            self.x = np.linspace(2200., altitude[0], num=15, endpoint=True)
        elif altitude == []:
            self.x = np.linspace(2200., 25000., num=15, endpoint=True)
        elif len(altitude) > 2 or altitude[0] > altitude[1]:
            print('bad altitude filter. Aborting')
            sys.exit()

        # Filtering by years or months or hours
        if years:
            self.dataframe = select_dataframe_by_year(self.dataframe, years)
            if months:
                self.dataframe = select_dataframe_by_month(self.dataframe, months)
            else:
                if epoch != 'all':
                    if self.observatory == 'north':
                        self.dataframe = select_new_epochs_dataframe_north(self.dataframe, epoch)
                    elif self.observatory == 'south':
                        self.dataframe = select_new_epochs_dataframe_south(self.dataframe, epoch)
                    else:
                        self.dataframe = select_new_epochs_dataframe_north(self.dataframe, epoch)
        elif months:
            self.dataframe = select_dataframe_by_month(self.dataframe, months)
        elif hours:
            self.dataframe = select_dataframe_by_hour(self.dataframe, hours)
        if not years and not months:
            if epoch != 'all':
                if self.observatory == 'north':
                    self.dataframe = select_new_epochs_dataframe_north(self.dataframe, epoch)
                elif self.observatory == 'south':
                    self.dataframe = select_new_epochs_dataframe_south(self.dataframe, epoch)
                else:
                    self.dataframe = select_new_epochs_dataframe_north(self.dataframe, epoch)

        if epoch_by_density:
            self.epoch = epoch_by_density_name
            selection = self.dataframe[(self.dataframe.n_exp > n_exp_minvalue) & (self.dataframe.n_exp < n_exp_maxvalue)
            & (self.dataframe.P == 125)]
            self.dataframe = self.dataframe[self.dataframe.MJD.isin(selection.MJD)]

        # Filtering by good weather conditions:
        if select_good_weather:
            self.dataframe['W'] = np.sqrt(self.dataframe.U ** 2. + self.dataframe.V ** 2.)
            h_cond = (self.dataframe.h > 2200.) & (self.dataframe.h < 2500.)
            gw_cond = (self.dataframe.RH < RH_lim) & (self.dataframe.W < W_lim * 1000. / 3600.)
            self.dataframe['datehour'] = self.dataframe.Date.astype(str) + self.dataframe.hour.astype(str)
            good_weather_dates = self.dataframe.datehour[(h_cond) & (gw_cond)]
            self.dataframe = self.dataframe[(self.dataframe.datehour.isin(good_weather_dates.tolist()))]

        # Various averaged values obtained after filtering (if no filter, averages are made over the whole dataframe)
        self.group_by_p = self.dataframe.groupby('P')
        self.h_avgs = avg_std_dataframe(self.group_by_p, 'h')
        self.n_avgs = avg_std_dataframe(self.group_by_p, 'n')
        self.Temp_avgs = avg_std_dataframe(self.group_by_p, 'Temp')
        self.wind_speed_avgs = avg_std_dataframe(self.group_by_p, 'wind_speed')
        self.wind_direction_avgs = avg_std_dataframe(self.group_by_p, 'wind_direction')
        self.RH_avgs = avg_std_dataframe(self.group_by_p, 'RH')
        self.n_exp_avgs = avg_std_dataframe(self.group_by_p, 'n_exp')

    def _compute_mass_density(self, air='moist', interpolate=False):
        """
        Uses the functions DensityMoistAir, MolarFractionWaterVapor and Compressibility
        from the LIDAR_analysis module humidity.py
        input:
            air: (str, optional) must be 'moist' or 'dry'

        :return: density [kg/m^3], std_dev(density) [kg/m^3], peak2peak_minus, peak2peak_plus
        """
        C = 415 # CO2 average global concentration in ppm
        rho_s = 1.225 # kg/m^3 standard air mass density (sea level, 15ºC)
        rho = []

        pbar = tqdm(total=len(self.dataframe.P))
        for i in np.arange(len(self.dataframe.P)):
            pbar.update(1)
            if air == 'moist':
                Xw = MolarFractionWaterVapor(self.dataframe.P.iloc[i], self.dataframe.Temp.iloc[i],
                                                  self.dataframe.RH.iloc[i])
                Z = Compressibility(self.dataframe.P.iloc[i], self.dataframe.Temp.iloc[i], Xw)
                rho.append(DensityMoistAir(self.dataframe.P.iloc[i] * 100., self.dataframe.Temp.iloc[i], Z, Xw, C))

            elif air == 'dry':
                Z = Compressibility(self.dataframe.P.iloc[i], self.dataframe.Temp.iloc[i], 0.0)
                rho.append(DensityMoistAir(self.dataframe.P.iloc[i]*100., self.dataframe.Temp.iloc[i], Z, 0.0, C))
            else:
                print('Wrong air condition. It must be "moist" or "dry". Aborting!')
                sys.exit()
        pbar.close()

        self.dataframe['n_mass_'+air] = rho
        self.dataframe['nexp_mass_' + air] = self.dataframe['n_mass_'+air] / rho_s * np.exp(self.dataframe.h / self.Hs)


    def _compute_diff_wrt_model(self):

        diff_with_magic = []
        diff_with_prod3 = []

        interpolated_density = self._interpolate_param_to_h('n_exp', self.x)[0]

        self._get_prod3sim_data()
        x = np.linspace(1000., 25000., num=15, endpoint=True)

        print("Computing the differences of the values of density:")
        for i in tqdm(np.arange(len(interpolated_density))):
            diff_with_magic.append((interpolated_density[i] - self.func_magic(x))
                                         / interpolated_density[i])
            diff_with_prod3.append((interpolated_density[i] - self.func_prod3(x))
                                         / interpolated_density[i])

        self.diff_with_magic = np.asarray(diff_with_magic)
        self.diff_with_prod3 = np.asarray(diff_with_prod3)
        self.diff_MAGIC = compute_averages_std(self.diff_with_magic)
        self.diff_PROD3 = compute_averages_std(self.diff_with_prod3)

    # =======================================================================================================
    # Plotting functions:
    # =======================================================================================================

    def plot_moist_dry_comparison(self, min_humidity=0.):

        self._compute_mass_density(air='dry')
        self._compute_mass_density(air='moist')
        
        dfmin_rh = self.dataframe[self.dataframe.RH > min_humidity]
        rel_dif = (dfmin_rh.nexp_mass_dry - dfmin_rh.nexp_mass_moist) * 100. / dfmin_rh.nexp_mass_dry
        
        dfmin_rh['rel_dif_n_mass'] = rel_dif
        group_by_P = dfmin_rh.groupby('P')
        height = avg_std_dataframe(group_by_P, 'h')
        
        rhow, erhow, emadrhow, pprhow, pmrhow = avg_std_dataframe(group_by_P, 'nexp_mass_moist')
        rhod, erhod, emadrhod, pprhod, pmrhod = avg_std_dataframe(group_by_P, 'nexp_mass_dry')
        rel_dif = avg_std_dataframe(group_by_P, 'rel_dif_n_mass')

        e_color_dry = '#1f77b4'
        e_color_moist = '#ff7f0e'
        
        fig, axs = plt.subplots(2,1,sharex=True)
        plt.subplots_adjust(hspace=0)
        
        axs[0].errorbar(height[0], rhod, yerr=emadrhod, fmt=':', color=e_color_dry, elinewidth=3, label=None)
        
        ebd = axs[0].errorbar(height[0], rhod, yerr=[pmrhod, pprhod], fmt='o', color=e_color_dry, capsize=0.5,
                              mec=e_color_dry, ms=2., label='ECMWF $\\rho_d$ (dry air)')
        
        axs[0].errorbar(height[0], rhow, yerr=emadrhow, fmt=':', color=e_color_moist, elinewidth=3, label=None)
        
        ebw = axs[0].errorbar(height[0], rhow, yerr=[pmrhow, pprhow], fmt='o', color=e_color_moist, capsize=0.5,
                              mec=e_color_moist, ms=2., label='ECMWF $\\rho_w$ (moist air)')
        
        ebd[-1][0].set_linestyle(':')
        ebw[-1][0].set_linestyle(':')

        axs[0].axvline(2000., ls='dotted')
        axs[0].set_ylim(10, 120)
        axs[0].set_ylabel('$\\rho / \\rho_{\\rm s}$ * exp(h/H$_{\\rm s}$) [kg m$^{-3}$]')
        axs[0].legend(loc='best')
        axs[0].axes.tick_params(direction='in')
        
        axs[1].errorbar(height[0], rel_dif[0], yerr=rel_dif[2], color='#8c564b', ms=2., label=None)
        ebrd = axs[1].errorbar(height[0], rel_dif[0], yerr=[rel_dif[4], rel_dif[3]], capsize=0.5, ms=2.,
                               color='#8c564b', mec='#8c564b', mfc='#8c564b', label='ECMWF')
        ebrd[-1][0].set_linestyle(':')

        axs[1].axvline(2000., ls='dotted')
        axs[1].legend(loc='best')
        axs[1].axes.tick_params(direction='inout', top='on')
        axs[1].set_xlabel('height [m]')
        axs[1].set_yscale('log')
        axs[1].set_ylabel('rel. diff [\\%]')
        fig.savefig('dry_vs_moist_air_density'+self.data_server+'_'+self.epoch+'.png', bbox_inches='tight', dpi=300)
        fig.show()

    def plot_average_at_15km(self, fmt='pdf'):
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

        xspan = ax.get_xlim()[1] - ax.get_xlim()[0]
        yspan = ax.get_ylim()[1] - ax.get_ylim()[0]

        for y in np.unique(self.dataframe.year):
            mjd_start_year = date2mjd(y, 1, 1, 0)
            mjd_half_year = date2mjd(y, 7, 1, 0)
            year_plot = y

            if 1 in np.unique(self.dataframe.month):
                ax.axvline(mjd_start_year, color='0.7', linewidth=1., linestyle='dotted', zorder=100)
                ax.text(mjd_start_year - xspan * 0.018, ax.get_ylim()[0] + 0.095 * yspan, str(year_plot), rotation=90,
                        color='0.7')

            if 7 in np.unique(self.dataframe.month):
                ax.axvline(mjd_half_year, color='0.7', linewidth=1., linestyle='dotted', zorder=100)
                ax.text(mjd_half_year - xspan * 0.018, ax.get_ylim()[0] + 0.095 * yspan, str(year_plot + 0.5),
                        rotation=90, color='0.7')

        ax.axhline(np.mean(density_at_15km), color='#336699', linestyle='solid', lw=1., zorder=10)

        ax.legend(loc='best', numpoints=1)
        ax.set_xlabel('MJD')
        ax.set_ylabel('$n_{\\rm day}/N_{\\rm s} \\cdot e^{(h/H_{\\rm s})}$')
        ax.set_ylim(np.min(density_at_15km) * 0.98, np.max(density_at_15km) * 1.02)
        ax.set_title('Density over time at h = 15 km (for %s months)' % (self.epoch))

        fig.savefig(self.output_plot_name + '_at_15_km.'+ fmt, bbox_inches='tight')
        fig.savefig(self.output_plot_name + '_at_15_km.png', bbox_inches='tight', dpi=300)

    def plot_differences_wrt_magic(self):
        """
        Function that plots the difference between density models w.r.t the MAGIC Winter model
        :return:
        """

        self._compute_diff_wrt_model()

        print('plotting averaged data values for selected epoch')
        fig = plt.figure()
        ax = fig.add_subplot(111)

        eb2 = ax.errorbar(self.x + 175., self.diff_MAGIC[0], yerr=[self.diff_MAGIC[3],
                                                                         self.diff_MAGIC[2]], fmt='o', color='b',
                          capsize=0.5, mec='b', ms=1., label=self.data_server)
        eb2[-1][0].set_linestyle(':')
        ax.errorbar(self.x + 175., self.diff_MAGIC[0], yerr=self.diff_MAGIC[1], fmt=':', color='b',
                    elinewidth=3.1)

        ax.axvline(2000., ls='dotted')
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

    def plot_differences_wrt_model(self, epochs, model, fmt='pdf'):

        fig = plt.figure()
        ax = fig.add_subplot(111)

        for e in epochs:
            self.get_data(e)

            self._compute_diff_wrt_model()

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

            color = next(ax._get_lines.prop_cycler)['color']
            eb2 = ax.errorbar(self.x, diff, yerr=ediff_pp, fmt='o', capsize=0.5, ms=1., label=e, color=color)
            eb2[-1][0].set_linestyle(':')
            ax.errorbar(self.x, diff, yerr=ediff, fmt=':', elinewidth=3.1, color=color)

        ax.axvline(2000., ls='dotted')
        ax.set_title('Relative Difference %s %s w.r.t %s model' % (self.data_server, self.observatory, model))
        ax.set_xlabel('h a.s.l. [m]')
        ax.set_ylabel('Rel. Difference')
        ax.set_xlim(0., 25100.)
        ax.set_ylim(-0.11, 0.09)
        ax.xaxis.set_minor_locator(MultipleLocator(1000))
        ax.xaxis.set_major_locator(MultipleLocator(2000))
        ax.yaxis.set_major_locator(MultipleLocator(0.02))
        ax.legend(loc='best', numpoints=1)
        ax.grid(which='both', axis='y', color='0.8')
        fig.savefig('differences_wrt_'+ model +'_' + self.output_plot_name +  '_' + self.observatory + '.' + fmt, bbox_inches='tight')
        fig.savefig('differences_wrt_'+ model +'_' + self.output_plot_name +  '_' + self.observatory + '.png', bbox_inches='tight', dpi=300)

    def plot_models_comparison(self, model=None, interpolate=False, fmt='pdf'):
        fig, ax = plt.subplots(2, 1, sharex=True)
        plt.subplots_adjust(hspace=0)
        if interpolate:
            raw_n_exp, avg_n_exp, e_n_exp, pp_n_exp, pm_n_exp = self._interpolate_param_to_h('n_exp', self.x)
            eb2 = ax[0].errorbar(self.x, avg_n_exp, yerr=[pm_n_exp, pp_n_exp],
                          fmt='o', color='b', capsize=0.5, mec='b', ms=1., label=self.data_server)
            eb2[-1][0].set_linestyle(':')
            ax[0].errorbar(self.x, avg_n_exp, yerr=e_n_exp, fmt=':', color='b', elinewidth=3., label=None)
            ax[1].plot(self.x, e_n_exp/avg_n_exp, 'o:', color='b', ms=2.)

        else:
            eb2 = ax[0].errorbar(self.h_avgs[0], self.n_exp_avgs[0], xerr=self.h_avgs[2], yerr=[self.n_exp_avgs[4],
                                                                                             self.n_exp_avgs[3]],
                              fmt='o', color='b', capsize=0.5, mec='b', ms=1., label=self.data_server)
            eb2[-1][0].set_linestyle(':')
            ax[0].errorbar(self.h_avgs[0], self.n_exp_avgs[0], xerr=self.h_avgs[2], yerr=self.n_exp_avgs[2], fmt=':',
                        color='b', elinewidth=3., label=None)
            ax[1].plot(self.h_avgs[0], self.n_exp_avgs[2]/self.n_exp_avgs[0], 'o:', color='b', ms=2.)

        if model == 'MW':
            ax[0].plot(self.heightmw * 1000., self.n_mw * np.exp(self.heightmw * 1000. / self.Hs), '-', color='grey',
                    label='MAGIC W')
        elif model == 'PROD3':
            self._get_prod3sim_data()
            ax[0].plot(self.hprod3 * 1000., self.n_prod3 * np.exp(self.hprod3 * 1000. / self.Hs), '-', color='0.2',
                    label='Prod3 ' + self.observatory)
        elif model == 'both':
            ax[0].plot(self.heightmw * 1000., self.n_mw * np.exp(self.heightmw * 1000. / self.Hs), '-', color='grey',
                    label='MAGIC W')
            self._get_prod3sim_data()
            ax[0].plot(self.hprod3 * 1000., self.n_prod3 * np.exp(self.hprod3 * 1000. / self.Hs), '-', color='0.2',
                    label='Prod3 ' + self.observatory)
        else:
            print('Wrong model. It must be "MW" for MAGIC Winter model, or "PROD3" for Paranal model or "both". '
                  '\n Exiting!')
            sys.exit()

        ax[0].axvline(2000., ls='dotted')
        ax[0].set_title(self.data_server + ' ' + self.observatory)
        ax[0].set_ylabel('$n_{\\rm day}/N_{\\rm s} \\cdot e^{(h/H_{\\rm s})}$')
        ax[0].set_xlim(0., 25100.)
        ax[0].set_ylim(0.4, 1.2)
        ax[0].xaxis.set_minor_locator(MultipleLocator(1000))
        ax[0].xaxis.set_major_locator(MultipleLocator(2000))
        ax[0].yaxis.set_major_locator(MultipleLocator(0.1))
        ax[0].legend(loc='best', numpoints=1)
        ax[0].grid(which='both', axis='y', color='0.8')
        ax[0].axes.tick_params(direction='in')

        ax[1].axvline(2000., ls='dotted')
        ax[1].axes.tick_params(direction='inout', top='on')
        ax[1].set_xlabel('h a.s.l. [m]')
        ax[1].grid(which='both', axis='y', color='0.8')
        ax[1].set_ylabel('std/$\\langle n_{\\rm day}/N_{\\rm s} \\cdot e^{(h/H_{\\rm s})} \\rangle$')

        fig.savefig('comparison_' + self.output_plot_name + '_' + self.observatory + '.' + fmt, bbox_inches='tight')
        fig.savefig('model_comparison_' + self.output_plot_name + '_' + self.observatory + '.png', bbox_inches='tight', dpi=300)

    def plot_epoch_comparison(self, param, epochs, hour=None, interpolate=False, plot_MW=False, plot_PROD3=False, format='png'):
        fig, ax = plt.subplots(2, 1, sharex=True)
        plt.subplots_adjust(hspace=0)
        for e in epochs:
            if hour:
                self.get_data(e, hours=[hour])
            else:
                self.get_data(e)
            color = next(ax[0]._get_lines.prop_cycler)['color']
            if interpolate:
                raw_n_exp, avg_n_exp, e_n_exp, pp_n_exp, pm_n_exp = self._interpolate_param_to_h('n_exp', self.x)
                eb2 = ax[0].errorbar(self.x, avg_n_exp, yerr=[pm_n_exp, pp_n_exp],
                              fmt='o', mec=None, color=color, capsize=0.5, ms=1., label=self.data_server+'\_'+e)
                eb2[-1][0].set_linestyle(':')
                ax[0].errorbar(self.x, avg_n_exp, yerr=e_n_exp, fmt=':', color=color, elinewidth=3., label=None)
                ax[1].plot(self.x, e_n_exp/avg_n_exp, 'o:', ms=2., color=color, label=self.data_server+'\_'+e)

            else:
                eb2 = ax[0].errorbar(self.h_avgs[0], self.n_exp_avgs[0], xerr=self.h_avgs[2], yerr=[self.n_exp_avgs[4],
                                                                                                    self.n_exp_avgs[3]],
                                     fmt='o', color=color, capsize=0.5, ms=1., label=self.data_server+'\_'+e)
                eb2[-1][0].set_linestyle(':')
                ax[0].errorbar(self.h_avgs[0], self.n_exp_avgs[0], xerr=self.h_avgs[2], yerr=self.n_exp_avgs[2],
                               fmt=':', color = color, elinewidth=3., label=None)
                ax[1].plot(self.h_avgs[0], self.n_exp_avgs[2]/self.n_exp_avgs[0], 'o:', ms=2., color=color,
                           label=self.data_server+'\_'+e)

        if plot_MW:
            ax[0].plot(self.heightmw * 1000., self.n_mw * np.exp(self.heightmw * 1000. / self.Hs), '-', color='grey',
                    label='MAGIC W')
        if plot_PROD3:
            self._get_prod3sim_data()
            ax[0].plot(self.hprod3 * 1000., self.n_prod3 * np.exp(self.hprod3 * 1000. / self.Hs), '-.', color='0.2',
                       label='Prod3 ' + self.observatory)

        ax[0].axvline(2000., ls='dotted')
        ax[0].set_title(self.data_server + self.observatory)
        ax[0].set_ylabel('$n_{\\rm day}/N_{\\rm s} \\cdot e^{(h/H_{\\rm s})}$')
        ax[0].set_xlim(0., 25100.)
        ax[0].set_ylim(0.4, 1.2)
        ax[0].xaxis.set_minor_locator(MultipleLocator(1000))
        ax[0].xaxis.set_major_locator(MultipleLocator(2000))
        ax[0].yaxis.set_major_locator(MultipleLocator(0.1))
        ax[0].legend(loc='best', numpoints=1)
        ax[0].grid(which='both', axis='y', color='0.8')
        ax[0].axes.tick_params(direction='in')

        ax[1].axvline(2000., ls='dotted')
        ax[1].axes.tick_params(direction='inout', top='on')
        ax[1].set_xlabel('h a.s.l. [m]')
        ax[1].set_ylabel('std/$\\langle n_{\\rm day}/N_{\\rm s} \\cdot e^{(h/H_{\\rm s})} \\rangle$')
        ax[1].set_ylim(0., 0.0449)
        ax[1].legend(loc='best', numpoints=1, ncol=2)
        ax[1].grid(which='both', axis='y', color='0.8')
        if hour:
            fig.savefig('epoch_comparison_' + self.output_plot_name + '_' + self.observatory + '_h' + str(hour) + '.' + format, bbox_inches='tight', dpi=300)
        else:
            fig.savefig('epoch_comparison_' + self.output_plot_name + '_' + self.observatory + '.' + format, bbox_inches='tight', dpi=300)

    # =======================================================================================================
    # printing functions:
    # =======================================================================================================

    def print_to_text_file(self):
        textfile = open(self.output_plot_name + '_to_text_file.txt', 'w')
        print('# x[i], avg, rms, p2p_m, p2p_p', file=textfile)
        for i in np.arange(len(self.x)):
            print(self.x[i], self.averages[0][i], self.averages[1][i], self.averages[3][i],
                  self.averages[2][i], self.diff_MAGIC[0][i], self.diff_MAGIC[1][i],
                  self.diff_MAGIC[3][i], self.diff_MAGIC[2][i], file=textfile)
        textfile.close()

    def _refractive_index(self, P, T, RH, wavelength):
        """Wrapper for Rayleigh.calculate_n()."""
        rayleigh = Rayleigh(wavelength, P, T, RH)
        return rayleigh.calculate_n()

    def write_corsika(self, outfile):
        """
        Write an output file in the style of a CORSIKA atmospheric configuration file:

        alt (km)     rho (g/cm^3)   thick (g/cm^2)   (n-1)
        """
        mbar2gcm2 = 1.019716213  # conversion from mbar (pressure in SI) to g/cm^2 (pressure in cgs)
        # Loschmidt constant: number density of particles in an ideal gas at STP in m^-3
        N0 = 2.079153e25  # Na divided by molar mass of air: 0.0289644
        height = np.array([1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20.,
                           21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 32., 34., 36., 38., 40., 42., 44., 46.,
                           48.])

        with open(outfile, 'w') as f:

            f.write("# Atmospheric Model ECMWF year/month/day   hour h\n")#.format(**datedict))
            f.write("#Col. #1          #2           #3            #4        [ #5 ]        [ #6 ]       [ # 7 ]\n")
            f.write("# Alt [km]    rho [g/cm^3] thick [g/cm^2]    n-1        T [k]       p [mbar]      pw / p\n")

            density = self.n_avgs[0].sort_index(ascending=False).values
            density /= N0*1.0E-3
            P = np.unique(self.dataframe.P)[::-1]

            #height = self.h_avgs[0].sort_index(ascending=False).values / 1.e3
            #Temp = self.Temp_avgs[0].sort_index(ascending=False).values
            #RH = self.RH_avgs[0].sort_index(ascending=False).values

            T0 = float(self._interpolate_simple(self.h_avgs[0].sort_index(ascending=False).values / 1e3,
                                                self.Temp_avgs[0].sort_index(ascending=False).values, 0.))
            RH0 = float(self._interpolate_simple(self.h_avgs[0].sort_index(ascending=False).values / 1e3,
                                                self.RH_avgs[0].sort_index(ascending=False).values, 0.))
            P0 = float(self._interpolate_simple(self.h_avgs[0].sort_index(ascending=False).values / 1e3, P, 0.))
            density0 = float(self._interpolate_simple(self.h_avgs[0].sort_index(ascending=False).values / 1e3, density,
                                                      0.))
            thick0 = P0 * mbar2gcm2

            pw0 = PartialPressureWaterVapor(T0, RH0) / P0

            Temp = self._interpolate_simple(self.h_avgs[0].sort_index(ascending=False).values / 1.e3,
                                                  self.Temp_avgs[0].sort_index(ascending=False).values, height)
            RH = self._interpolate_simple(self.h_avgs[0].sort_index(ascending=False).values / 1.e3,
                                                self.RH_avgs[0].sort_index(ascending=False).values, height)
            RH[RH < 0.] = 1.e-4
            Pressure = self._interpolate_simple(self.h_avgs[0].sort_index(ascending=False).values / 1.e3, P, height)
            density = self._interpolate_simple(self.h_avgs[0].sort_index(ascending=False).values / 1.e3, density,
                                                     height)
            thick = Pressure * mbar2gcm2
            pwp = PartialPressureWaterVapor(Temp, RH) / Pressure

            nm0 = self._refractive_index(P0, T0, RH0, 350.) - 1.
            outdict = {'height': 0.000, 'rho': density0, 'thick': thick0, 'nm1': nm0, 'T':T0, 'p': P0,
                       'pw/p': pw0}
            f.write("  {height:7.3f}     {rho:5.5E}  {thick:5.5E}  {nm1:5.5E}  {T:5.5E}  {p:5.5E}  {pw/p:5.5E}\n".format(**outdict))

            for i in np.arange(len(height)):
                nm1 = self._refractive_index(Pressure[i], Temp[i], RH[i], 350) - 1

                outdict = {'height': height[i], 'rho': density[i], 'thick': thick[i], 'nm1': nm1, 'T': Temp[i],
                           'p': Pressure[i], 'pw/p': pwp[i]}
                f.write("  {height:7.3f}     {rho:5.5E}  {thick:5.5E}  {nm1:5.5E}  {T:5.5E}  {p:5.5E}  {pw/p:5.5E}\n".format(**outdict))

            # concatenate the dummy values of MagicWinter starting from 50.0 km
            f.write("   50.000     1.09738e-06  8.72656e-01   2.53289e-07   266.34   8.38955e-01  0.00000e+00\n")
            f.write("   55.000     5.99974e-07  4.61036e-01   1.38481e-07   257.19   4.42930e-01  0.00000e+00\n")
            f.write("   60.000     3.25544e-07  2.36175e-01   7.51395e-08   242.81   2.26896e-01  0.00000e+00\n")
            f.write("   65.000     1.70152e-07  1.15918e-01   3.92732e-08   227.93   1.11324e-01  0.00000e+00\n")
            f.write("   70.000     8.43368e-08  5.45084e-02   1.94660e-08   215.90   5.22651e-02  0.00000e+00\n")
            f.write("   75.000     3.95973e-08  2.48012e-02   9.13953e-09   208.66   2.37169e-02  0.00000e+00\n")
            f.write("   80.000     1.79635e-08  1.10899e-02   4.14618e-09   205.11   1.05760e-02  0.00000e+00\n")
            f.write("   85.000     8.03691e-09  4.91583e-03   1.85502e-09   202.12   4.66284e-03  0.00000e+00\n")
            f.write("   90.000     3.59602e-09  2.15599e-03   8.30003e-10   196.26   2.02583e-03  0.00000e+00\n")
            f.write("   95.000     1.59871e-09  9.21029e-04   3.69000e-10   187.55   8.60656e-04  0.00000e+00\n")
            f.write("  100.000     6.73608e-10  3.82814e-04   1.55477e-10   185.38   3.58448e-04  0.00000e+00\n")
            f.write("  105.000     2.69097e-10  1.61973e-04   6.21108e-11   197.19   1.52311e-04  0.00000e+00\n")
            f.write("  110.000     1.09021e-10  7.37110e-05   2.51634e-11   224.14   7.01416e-05  0.00000e+00\n")
            f.write("  115.000     4.71300e-11  3.70559e-05   1.08782e-11   268.51   3.63251e-05  0.00000e+00\n")
            f.write("  120.000     2.23479e-11  2.05900e-05   5.15817e-12   333.43   2.13890e-05  0.00000e+00\n")
