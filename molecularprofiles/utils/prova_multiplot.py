from molecularprofiles.molecularprofiles import *
from LIDAR_Analysis.humidity import *
epoch='all'
gdas = MolecularProfile('/home/pmunar/feina/software/pandas_branch/molecularprofiles/molecularprofiles/test_data/gdas_merged_2012-01_2017-08_north.txt',data_server='GDAS', observatory='south')
gdas.get_data(epoch=epoch)
ecmwf = MolecularProfile('/home/pmunar/feina/software/pandas_branch/molecularprofiles/molecularprofiles/test_data/ecmwf_interim_merged_2012-01_2017-12_north.txt', data_server='ECMWF', observatory='north')
ecmwf.get_data(epoch=epoch)


fig, ax = plt.subplots(nrows=3, ncols=4, figsize=(6, 6*(np.sqrt(5.)- 1.0)/2.), sharey=True, sharex=True)
years = np.unique(ecmwf.dataframe.year)
years = years.astype(int).tolist()
years = [2012,2013,2014,2015,2016]
months = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
          'November', 'December']
for y in years:
    m = 1
    for a in ax.flatten():
        print(y, m)
        ecmwf.get_data(years=[y], months=[m], select_good_weather=True, RH_lim=95., W_lim=50.)
        a.plot(ecmwf.h_avgs[0], ecmwf.n_exp_avgs[2]/ecmwf.n_exp_avgs[0], label=str(y))
        a.axvline(2200., ls=':', color='k')
        a.set_title(months[m-1])

        m += 1


ecmwf.get_data('winter')

for c in range(4):
    ax[0,c].plot(ecmwf.h_avgs[0], ecmwf.n_exp_avgs[2] / ecmwf.n_exp_avgs[0], label='winter', lw=2.)
ax[2,3].plot(ecmwf.h_avgs[0], ecmwf.n_exp_avgs[2] / ecmwf.n_exp_avgs[0], label='winter', lw=2.)


ecmwf.get_data('intermediate')
ax[1,0].plot(ecmwf.h_avgs[0], ecmwf.n_exp_avgs[2] / ecmwf.n_exp_avgs[0], label='intermediate', lw=2.)
ax[1,1].plot(ecmwf.h_avgs[0], ecmwf.n_exp_avgs[2] / ecmwf.n_exp_avgs[0], label='intermediate', lw=2.)
ax[2,1].plot(ecmwf.h_avgs[0], ecmwf.n_exp_avgs[2] / ecmwf.n_exp_avgs[0], label='intermediate', lw=2.)
ax[2,2].plot(ecmwf.h_avgs[0], ecmwf.n_exp_avgs[2] / ecmwf.n_exp_avgs[0], label='intermediate', lw=2.)

ecmwf.get_data('summer')
ax[1,2].plot(ecmwf.h_avgs[0], ecmwf.n_exp_avgs[2] / ecmwf.n_exp_avgs[0], label='summer', lw=2.)
ax[1,3].plot(ecmwf.h_avgs[0], ecmwf.n_exp_avgs[2] / ecmwf.n_exp_avgs[0], label='summer', lw=2.)
ax[2,0].plot(ecmwf.h_avgs[0], ecmwf.n_exp_avgs[2] / ecmwf.n_exp_avgs[0], label='summer', lw=2.)
for a in ax.flatten():
#     a.plot(ecmwf.h_avgs[0], ecmwf.n_exp_avgs[2] / ecmwf.n_exp_avgs[0], label='all', lw=2.)
    a.legend(ncol=2, loc='upper right')
    a.set_xlim(0., 25000.)
for c in range(4):
    ax[2,c].set_xlabel('height [m]')
for c in range(3):
    ax[c,0].set_ylabel('mad(n\_exp)/mean(n\_exp)')

fig.show()

plt.hlines(6.0, mean_lh30 - sigma3_lh30, mean_lh30 + sigma3_lh30)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6*(np.sqrt(5.)- 1.0)/2.), sharey=True)

epochs = ['winter', 'summer', 'intermediate']

ax.set_color_cycle(clist=['#1f77b4', '#ff7f0e', '#2ca02c'])

#for y in years:
for e in epochs:
#    print('year %i, epoch %s' %(y, e))
    months = get_epoch(e)

    ecmwf.get_data(e, select_good_weather=True, RH_lim=95., W_lim=50.)
    density_at_15km = ecmwf._interpolate_param_to_h('n_exp', 15000.)
    ax.plot(np.unique(ecmwf.dataframe.MJD), density_at_15km, 'o', markersize=1.2,
            label=str(e), alpha=0.8)

ax.set_xlabel('MJD')
ax.set_ylabel('$n_{\\rm day}/N_{\\rm s} \\cdot e^{(h/H_{\\rm s})}$')
#        ax.set_title(y)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[0:3], labels[0:3], frameon=True)

xspan = ax.get_xlim()[1] - ax.get_xlim()[0]
yspan = ax.get_ylim()[1] - ax.get_ylim()[0]

for y in years:
    mjd_start_year = date2mjd(y, 1, 1, 0)
    mjd_half_year = date2mjd(y, 7, 1, 0)
    year_plot = y

    print('true 1')
    ax.axvline(mjd_start_year, color='0.7', linewidth=1., linestyle='dotted', zorder=100)
    ax.text(mjd_start_year - xspan * 0.018, ax.get_ylim()[0] + 0.095 * yspan, str(year_plot), rotation=90, color='0.3')

    print('true 2')
    ax.axvline(mjd_half_year, color='0.7', linewidth=1., linestyle='dotted', zorder=100)
    ax.text(mjd_half_year - xspan * 0.018, ax.get_ylim()[0] + 0.095 * yspan, str(year_plot + 0.5), rotation=90,
            color='0.3')

fig.show()

file = open('list_of_2016txts.txt', 'r')

line = file.readline()[:-1]
epoch = 'all'
while line:
    ecmwf = MolecularProfile(line, data_server='ECMWF', observatory='Krakow')
    ecmwf.get_data(epoch=epoch)
    ecmwf.write_corsika(line.split('.')[0]+'_corsika_input.txt')
    line = file.readline()[:-1]


