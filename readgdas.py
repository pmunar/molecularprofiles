def read_gdas(file_gdas, epoch_text):
    import numpy as np
    import sys
    from grib_utils import mjd2date
    from grib_utils import get_epoch

    print("loading and selecting gdas data")
    gdas_file = open(file_gdas)
    mjd_gdas, h_gdas, nns_gdas, p_gdas = np.loadtxt(gdas_file, usecols=(0, 1, 2, 5), unpack=True)

    epoch = get_epoch(epoch_text)

    year_gdas = []
    month_gdas = []
    day_gdas = []
    hour_gdas = []
    for i in np.arange(len(mjd_gdas)):
        # Percentage counter bar
        sys.stdout.write('\r')
        # the exact output you're looking for:
        k = int(i * 101/len(mjd_gdas))
        sys.stdout.write("[%-100s] %d%%" % ('=' * k, k))
        sys.stdout.flush()
        # -------------------------------
        date_gdas = mjd2date(mjd_gdas[i])
        year_gdas.append(date_gdas[0])
        month_gdas.append(date_gdas[1])
        day_gdas.append(date_gdas[2])
        hour_gdas.append(date_gdas[3])
    print('\n')
    year_gdas = np.asarray(year_gdas)
    month_gdas = np.asarray(month_gdas)
    day_gdas = np.asarray(day_gdas)
    hour_gdas = np.asarray(hour_gdas)

    condition_1 = (mjd_gdas % 1 != 0.125) & (mjd_gdas % 1 != 0.375) & (mjd_gdas % 1 != 0.625) & \
                  (mjd_gdas % 1 != 0.875) & (p_gdas < 1001)

    h_gdas = h_gdas[condition_1] * 1000.
    nns_gdas = nns_gdas[condition_1]
    year_gdas = year_gdas[condition_1]
    month_gdas = month_gdas[condition_1]
    day_gdas = day_gdas[condition_1]
    hour_gdas = hour_gdas[condition_1]
    tmp_p_gdas = p_gdas[condition_1]
    mjd_gdas = mjd_gdas[condition_1]

    p_gdas = tmp_p_gdas
    # Selecting now the data corresponding to the "epoch" only:


    if epoch_text == 'winter':
        mjd_gdas = mjd_gdas[(month_gdas == epoch[0]) | (month_gdas == epoch[1]) | (month_gdas == epoch[2]) |
                            (month_gdas == epoch[3]) | (month_gdas == epoch[4])]
        h_gdas = h_gdas[(month_gdas == epoch[0]) | (month_gdas == epoch[1]) | (month_gdas == epoch[2]) |
                        (month_gdas == epoch[3]) | (month_gdas == epoch[4])]
        nns_gdas = nns_gdas[(month_gdas == epoch[0]) | (month_gdas == epoch[1]) | (month_gdas == epoch[2]) |
                            (month_gdas == epoch[3]) | (month_gdas == epoch[4])]
        p_gdas = p_gdas[(month_gdas == epoch[0]) | (month_gdas == epoch[1]) | (month_gdas == epoch[2]) |
                        (month_gdas == epoch[3]) | (month_gdas == epoch[4])]

    if epoch_text == 'summer':
        mjd_gdas = mjd_gdas[(month_gdas == epoch[0]) | (month_gdas == epoch[1]) | (month_gdas == epoch[2])]
        h_gdas = h_gdas[(month_gdas == epoch[0]) | (month_gdas == epoch[1]) | (month_gdas == epoch[2])]
        nns_gdas = nns_gdas[(month_gdas == epoch[0]) | (month_gdas == epoch[1]) | (month_gdas == epoch[2])]
        p_gdas = p_gdas[(month_gdas == epoch[0]) | (month_gdas == epoch[1]) | (month_gdas == epoch[2])]

    if epoch_text == 'intermediate':
        mjd_gdas = mjd_gdas[(month_gdas == epoch[0]) | (month_gdas == epoch[1]) | (month_gdas == epoch[2]) |
                            (month_gdas == epoch[3])]
        h_gdas = h_gdas[(month_gdas == epoch[0]) | (month_gdas == epoch[1]) | (month_gdas == epoch[2]) |
                        (month_gdas == epoch[3])]
        nns_gdas = nns_gdas[(month_gdas == epoch[0]) | (month_gdas == epoch[1]) | (month_gdas == epoch[2]) |
                            (month_gdas == epoch[3])]
        p_gdas = p_gdas[(month_gdas == epoch[0]) | (month_gdas == epoch[1]) | (month_gdas == epoch[2]) |
                        (month_gdas == epoch[3])]

    if epoch_text == 'all':
        mjd_gdas = mjd_gdas
        h_gdas = h_gdas
        nns_gdas = nns_gdas
        p_gdas = p_gdas

    return mjd_gdas, year_gdas, month_gdas, day_gdas, hour_gdas, h_gdas, nns_gdas, p_gdas


def create_mjd_list(mjd_ecmwf):
    import numpy as np
    # This is to create a list of the observed dates (mjd's in this case). This is because mjd_ecmwf is in a format
    # shuch that there is a row for every preassure level, and hence there are as many rows as pressure levels with the
    # same MDJ. We want every different MJD, and this function does that. Reads mjd_ecmwf and records only those
    # elements that are different from each other.

    mjd_list = []
    mjd_test = 0
    for i in np.arange(len(mjd_ecmwf)):
        if mjd_test != mjd_ecmwf[i]:
            mjd_list.append(mjd_ecmwf[i])
            mjd_test = mjd_ecmwf[i]
    return np.asarray(mjd_list)


def rearange_gdas(file_gdas, epoch_text, mjd_ecmwf):
    import numpy as np
    from grib_utils import mjd2date, get_epoch

    mjd_gdas, year_gdas, month_gdas, day_gdas, hour_gdas, h_gdas, nns_gdas, p_gdas = read_gdas(file_gdas, epoch_text)

    # We have to invert the data in each day, without inverting the order of the days,
    # bacause the values in GDAS data are in inverse order with respect to ECMWF data:
    # Without this, we would be comparing the Plevel 1 from GDAS with the Plevel 1000
    # from ECMWF for each day.

    mjd_list = create_mjd_list(mjd_ecmwf)

    n_mjd_gdas = []
    n_p_gdas = []
    n_nns_gdas = []
    n_h_gdas = []
    n_year_gdas = []
    n_month_gdas = []
    n_day_gdas = []
    n_hour_gdas = []
    for d in mjd_list:
        mjds = mjd_gdas[mjd_gdas == d].tolist()[::-1]
        ps = p_gdas[mjd_gdas == d].tolist()[::-1]
        hs = h_gdas[mjd_gdas == d].tolist()[::-1]
        nnss = nns_gdas[mjd_gdas == d].tolist()[::-1]
        years = year_gdas[mjd_gdas == d].tolist()[::-1]
        months = month_gdas[mjd_gdas == d].tolist()[::-1]
        days = day_gdas[mjd_gdas == d].tolist()[::-1]
        hours = hour_gdas[mjd_gdas == d].tolist()[::-1]

        n_mjd_gdas.extend(mjds)
        n_p_gdas.extend(ps)
        n_h_gdas.extend(hs)
        n_nns_gdas.extend(nnss)
        n_year_gdas.extend(years)
        n_month_gdas.extend(months)
        n_day_gdas.extend(days)
        n_hour_gdas.extend(hours)

    h_gdas = np.asarray(n_h_gdas)
    p_gdas = np.asarray(n_p_gdas)
    nns_gdas = np.asarray(n_nns_gdas)
    mjd_gdas = np.asarray(n_mjd_gdas)
    year_gdas = np.asarray(n_year_gdas)
    month_gdas = np.asarray(n_month_gdas)
    day_gdas = np.asarray(n_day_gdas)
    hour_gdas = np.asarray(n_hour_gdas)
    return mjd_gdas, year_gdas, month_gdas, day_gdas, hour_gdas, h_gdas, nns_gdas, p_gdas
