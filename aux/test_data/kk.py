mwd = ecmwf.dataframe.wind_direction
mwd[mwd < 0.] += 360.
wd = 270 - mwd
wd[wd < 0.] += 360.


fig = plt.figure()
for h in [2000., 6000., 10000., 15000., 20000., 25000.]:
    wind_speed_at_h = ecmwf._interpolate_param_to_h('wind_speed', h)
    wind_dir_at_h = ecmwf._interpolate_param_to_h('wind_direction', h)
    ax = plt.subplot(111)

    cb = ax.scatter(np.unique(ecmwf.dataframe.MJD), wind_dir_at_h, c = wind_speed_at_h, marker='o', s=1.5, alpha=0.4)
    ccb = plt.colorbar(cb, ticks=np.arange(0,100,10).tolist(), boundaries=np.arange(0.,101., 0.5))
    ccb.set_clim(0., 100.)
    ccb.set_label('wind speed [m/s]')

    ax.set_xlabel('MJD')
    ax.set_ylabel('wind direction [$^\\circ$]')
    ax.set_title('winds at h = %s' %(h))
    ax.set_ylim(-60, 400.)

    xspan = ax.get_xlim()[1] - ax.get_xlim()[0]
    yspan = ax.get_ylim()[1] - ax.get_ylim()[0]

    for y in np.unique(ecmwf.dataframe.year):
        mjd_start_year = date2mjd(y, 1, 1, 0)
        mjd_half_year = date2mjd(y, 7, 1, 0)
        year_plot = y
        if 1 in np.unique(ecmwf.dataframe.month):
            ax.axvline(mjd_start_year, color='0.7', linewidth=1., linestyle='dotted', zorder=100)
            ax.text(mjd_start_year - xspan * 0.018, ax.get_ylim()[0] + 0.095 * yspan, str(year_plot), rotation=90,
                    color='0.3')
        if 7 in np.unique(ecmwf.dataframe.month):
            ax.axvline(mjd_half_year, color='0.7', linewidth=1., linestyle='dotted', zorder=100)
            ax.text(mjd_half_year - xspan * 0.018, ax.get_ylim()[0] + 0.095 * yspan, str(year_plot + 0.5), rotation=90,
                    color='0.3')




    fig.savefig('wind_dir_at_h'+str(h)+'.png', dpi=300, bbox_inches='tight')
    fig.clf()


fig = plt.figure()
for h in [2000., 6000., 10000., 15000., 20000., 25000.]:
    wind_speed_at_h = ecmwf._interpolate_param_to_h('wind_speed', h)
    wind_dir_at_h = ecmwf._interpolate_param_to_h('wind_direction', h)
    ax = plt.subplot(111, projection='polar')
    cb = ax.scatter(wind_dir_at_h*np.pi/180., wind_speed_at_h, c=ecmwf.dataframe.month[ecmwf.dataframe.P == 1.], marker='o', s=1.5, alpha=0.4)
    ccb = plt.colorbar(cb)
    ccb.set_label('month of year')
    ax.set_rmax(100.)
    ax.set_rlabel_position(-22.5)
    ax.grid(True)
    fig.savefig('polar_plot_wind_h'+str(h)+'.png', dpi=300, bbox_inches='tight')
    fig.clf()

epochs = ['winter', 'summer', 'intermediate']
for h in np.arange(2000., 20000., 500).tolist():
    fig = plt.figure()
    ax = plt.subplot(111, projection='polar')
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location("N")
    for e in epochs:
        ecmwf.get_data(epoch=e)
        wind_speed_at_h = ecmwf._interpolate_param_to_h('wind_speed', h)
        wind_dir_at_h = ecmwf._interpolate_param_to_h('wind_direction', h)
        ax.scatter(wind_dir_at_h*np.pi/180., wind_speed_at_h,
                   marker='o', s=1.5, alpha=0.4, label=e, color = next(ax._get_lines.prop_cycler)['color'])
        #ccb = plt.colorbar(cb)
        #ccb.set_label('month of year')
    ax.set_rmax(100.)
    ax.set_rlabel_position(-22.5)
    ax.set_title('Wind direction and speed at altitude = %s'%(h))
    ax.legend(frameon=True, fancybox=True, framealpha=0.7)
    ax.grid(True)
    fig.savefig('polar_plot_wind_h'+str(h)+'.pdf', bbox_inches='tight')
#   fig.show()
    fig.clf()


def create_wind_speed_dataframe(df):
    wd_centre_bins = np.arange(7.5,360, 15)
    ws_hist = []
    for d in wd_centre_bins:
        ws_hist.append(np.histogram(df.wind_speed[(df.wind_direction >= d-7.5) & (df.wind_direction < d + 7.5)],
                                    bins=[0,5,10,20,30,40,50,100])[0])

    df_winds = pd.DataFrame(columns=['direction', '0-5', '5-10', '10-20', '20-30', '30-40', '40-50', '> 50'])
    ws_new_list = []
    for j in range(len(ws_hist[0])):
        li = []
        for i in range(len(ws_hist)):
            li.append(ws_hist[i][j])
        ws_new_list.append(li)

    for i,j in zip(df_winds.keys()[1:], range(len(ws_new_list))):
        df_winds[i] = ws_new_list[j]
    df_winds_normalized = df_winds.div(df_winds.sum(axis=1), axis=0)
    df_winds_normalized['direction'] = wd_centre_bins
    return df_winds_normalized

def plot_wind_rose(df, name_tag='my_wind_rose'):
    data = []
    counter = 0
    for col in df.columns:
        if col != 'direction':
            data.append(go.Area(t=df['direction'], r=df[col],
                                marker=dict(color=cl.scales['9']['seq']['YlGnBu'][counter]), name=col+' m/s'))
            counter+=1
    print(data)

    fig = go.Figure(data=data, layout=go.Layout(orientation=270., barmode='stack'))
    pio.write_image(fig, name_tag+'_wind_speed_rose.pdf')
    plotly.offline.plot(fig)
