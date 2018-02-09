from molecularprofiles.molecularprofiles import *

epoch='intermediate'

gdas = MolecularProfile('/home/pmunar/feina/meteo_data/grib_files/gdas/all_2012-2017_merged_file.txt',
                        data_server='GDAS', )
gdas.get_data(epoch=epoch)
ecmwf = MolecularProfile('/home/pmunar/feina/meteo_data/grib_files/ecmwf/txt_data_files/interim_2012_01-2017_08_merged_file.txt')
ecmwf.get_data(epoch=epoch)



ecmwf.compute_mass_density(air='dry')
ecmwf.compute_mass_density(air='moist')
gdas.compute_mass_density(air='moist')
gdas.compute_mass_density(air='dry')


edf80 = ecmwf.dataframe[ecmwf.dataframe.RH > 80.]
gdf80 = gdas.dataframe[gdas.dataframe.RH > 80.]

erel_dif = (edf80.nexp_mass_dry - edf80.nexp_mass_moist) * 100. / edf80.nexp_mass_dry
grel_dif = (gdf80.nexp_mass_dry - gdf80.nexp_mass_moist) * 100. / gdf80.nexp_mass_dry

edf80['rel_dif_n_mass'] = erel_dif
gdf80['rel_dif_n_mass'] = grel_dif

egp = edf80.groupby('P')
ggp = gdf80.groupby('P')

eh = avg_std_dataframe(egp, 'h')
gh = avg_std_dataframe(ggp, 'h')

erhow, eerhow, pperhow, pmerhow = avg_std_dataframe(egp, 'n_mass_moist')
erhod, eerhod, pperhod, pmerhod = avg_std_dataframe(egp, 'n_mass_dry')
grhod, egrhod, ppgrhod, pmgrhod = avg_std_dataframe(ggp, 'n_mass_dry')
grhow, egrhow, ppgrhow, pmgrhow = avg_std_dataframe(ggp, 'n_mass_moist')
erel_dif = avg_std_dataframe(egp, 'rel_dif_n_mass')
grel_dif = avg_std_dataframe(ggp, 'rel_dif_n_mass')



color_dry = '#1f77b4'
g_color_moist = '#ff7f0e'
e_color_dry = '#1f77b4'
e_color_moist = '#ff7f0e'
g_color_dry = '#2ca02c'
g_color_moist = '#d62728'

fig, axs = plt.subplots(2,1,sharex=True)
plt.subplots_adjust(hspace=0)
axs[0].errorbar(eh[0], erhod, yerr=eerhod, fmt=':', color=e_color_dry, elinewidth=3, label=None)

ebd = axs[0].errorbar(eh[0], erhod, yerr=[pmerhod, pperhod], fmt='o', color=e_color_dry, capsize=0.5, mec=e_color_dry, ms=2., label='ECMWF $\\rho_d$ (dry air)')

axs[0].errorbar(eh[0], erhow, yerr=eerhow, fmt=':', color=e_color_moist, elinewidth=3, label=None)

ebw = axs[0].errorbar(eh[0], erhow, yerr=[pmerhow, pperhow], fmt='o', color=e_color_moist, capsize=0.5, mec=e_color_moist, ms=2., label='ECMWF $\\rho_w$ (moist air)')

axs[0].errorbar(gh[0], grhow, yerr=egrhow, fmt=':', color=g_color_moist, elinewidth=3, label=None)

gbw = axs[0].errorbar(gh[0], grhow, yerr=[pmgrhow, ppgrhow], fmt='o', color=g_color_moist, capsize=0.5, mec=g_color_moist, ms=2., label='GDAS $\\rho_w$ (moist air)')

axs[0].errorbar(gh[0], grhod, yerr=egrhod, fmt=':', color=g_color_dry, elinewidth=3, label=None)
gbd = axs[0].errorbar(gh[0], grhod, yerr=[pmgrhod, ppgrhod], fmt='o', color=g_color_dry, capsize=0.5, mec=g_color_dry, ms=2., label='GDAS $\\rho_d$ (dry air)')
ebd[-1][0].set_linestyle(':')
ebw[-1][0].set_linestyle(':')
gbd[-1][0].set_linestyle(':')
gbw[-1][0].set_linestyle(':')

axs[0].set_ylabel('$\\rho$ * exp(h/H$_{\\rm s}$) [kg m$^{-3}$]')
axs[0].legend(loc='best')
axs[0].axes.tick_params(direction='in')

axs[1].errorbar(eh[0], erel_dif[0], yerr=erel_dif[1], color='#8c564b', ms=2., label=None)
ebrd = axs[1].errorbar(eh[0], erel_dif[0], yerr=[erel_dif[3], erel_dif[2]], capsize=0.5, ms=2., color='#8c564b', mec='#8c564b', mfc='#8c564b', label='ECMWF')
axs[1].errorbar(gh[0], grel_dif[0], yerr=grel_dif[1], color='#17becf', ms=2., label=None)
gbrd = axs[1].errorbar(gh[0], grel_dif[0], yerr=[grel_dif[3], grel_dif[2]], capsize=0.5, ms=2., color='#17becf', mec='#17becf', mfc='#17becf', label='GDAS')
ebrd[-1][0].set_linestyle(':')
gbrd[-1][0].set_linestyle(':')

axs[1].legend(loc='best')
axs[1].axes.tick_params(direction='inout', top='on')
axs[1].set_xlabel('height [m]')
axs[1].set_yscale('log')
axs[1].set_ylabel('rel. diff [\\%]')
fig.savefig('dry_vs_moist_air_density_'+epoch+'.png', bbox_inches='tight', dpi=300)
fig.show()


fig, axs = plt.subplots(2,1,sharex=True)
plt.subplots_adjust(hspace=0)

e_interpolated = ecmwf._interpolate_param_to_h('n_exp', ecmwf.x)
g_interpolated = gdas._interpolate_param_to_h('n_exp', gdas.x)

rel_dif = (e_interpolated[0] - g_interpolated[0]) / e_interpolated[0] * 100
av_rel_dif = compute_averages_std(rel_dif)

fig, axs = plt.subplots(2,1,sharex=True)
plt.subplots_adjust(hspace=0)
color_ecmwf = '#1f77b4'
color_gdas = '#ff7f0e'
axs[0].errorbar(ecmwf.x, e_interpolated[1], yerr=e_interpolated[2], ms=2., label=None, color=color_ecmwf, fmt=':',
                elinewidth=3)
eb1 = axs[0].errorbar(ecmwf.x, e_interpolated[1], yerr=[e_interpolated[4], e_interpolated[3]], capsize=0.5, ms=2.,
                      color=color_ecmwf, mec=color_ecmwf, mfc=color_ecmwf, label='ECMWF')

axs[0].errorbar(gdas.x, g_interpolated[1], yerr=g_interpolated[2], ms=2, label=None, color=color_gdas, fmt=':',
                elinewidth=3)
eb2 = axs[0].errorbar(gdas.x, g_interpolated[1], yerr=[g_interpolated[4], g_interpolated[3]], capsize=0.5, ms=2.,
                      color=color_gdas, mec=color_gdas, mfc=color_gdas, label='GDAS')

eb1[-1][0].set_linestyle(':')
eb2[-1][0].set_linestyle(':')

axs[1].errorbar(ecmwf.x, av_rel_dif[0], yerr=av_rel_dif[1], ms=2., label=None, color='#d62728', fmt=':', elinewidth=3)
eb3 = axs[1].errorbar(ecmwf.x, av_rel_dif[0], yerr=[av_rel_dif[3], av_rel_dif[2]], capsize=0.5, ms=2., color='#d62728',
                      mec='#d62728' , mfc='#d62728', label='ECMWF - GDAS')
eb3[-1][0].set_linestyle(':')

axs[1].legend(loc='best')
axs[1].axes.tick_params(direction='inout', top='on')
axs[1].set_xlabel('height [m]')
axs[1].set_ylabel('rel. diff [\\%]')
axs[0].set_ylabel('$\\rho$ * exp(h/H$_{\\rm s}$) [kg m$^{-3}$]')
axs[0].legend(loc='best')
axs[0].axes.tick_params(direction='in')
fig.savefig('ecmwf_vs_gdas_air_particle_density_'+epoch+'.png', bbox_inches='tight', dpi=300)
fig.show()
