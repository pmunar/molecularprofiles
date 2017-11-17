def read_ecmwf(file_ecmwf, epoch_text):
    import numpy as np
    from time import sleep
    import sys
    from grib_utils import date2mjd, get_epoch
    global months, month_ecmwf
    print("loading and selecting ecmwf data")
    ecmwf_file = open(file_ecmwf)
    date_ecmwf, year_ecmwf, month_ecmwf, day_ecmwf, hour_ecmwf, p_ecmwf, T_ecmwf, geop_ecmwf, h_ecmwf, n_ecmwf = \
        np.loadtxt(ecmwf_file, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), unpack=True)

    mjd_ecmwf = []

    for i in np.arange(len(date_ecmwf)):
        # Percentage counter bar
        sys.stdout.write('\r')
        k = int(i * 101/len(date_ecmwf))
        sys.stdout.write("[%-100s] %d%%" % ('=' * k, k))
        sys.stdout.flush()
        # ----------------------------
        mjd_ecmwf.append(date2mjd(year_ecmwf[i], month_ecmwf[i], day_ecmwf[i], hour_ecmwf[i]))
    print('\n')
    mjd_ecmwf = np.asarray(mjd_ecmwf)
    epoch = get_epoch(epoch_text)

    # mjd_ecmwf = mjd_ecmwf[month_ecmwf == epoch]
    # year_ecmwf = year_ecmwf[month_ecmwf == epoch]
    # day_ecmwf = date_ecmwf[month_ecmwf == epoch]
    # hour_ecmwf = hour_ecmwf[month_ecmwf == epoch]
    # h_ecmwf   = h_ecmwf[month_ecmwf == epoch]
    # dn_ecmwf   = n_ecmwf[month_ecmwf == epoch]
    # p_ecmwf   = p_ecmwf[month_ecmwf == epoch]
    # month_ecmwf = month_ecmwf[month_ecmwf == epoch]

    if epoch_text == 'winter':
        mjd_ecmwf = mjd_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                              (month_ecmwf == epoch[3]) | (month_ecmwf == epoch[4])]
        h_ecmwf = h_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                          (month_ecmwf == epoch[3]) | (month_ecmwf == epoch[4])]
        n_ecmwf = n_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                          (month_ecmwf == epoch[3]) | (month_ecmwf == epoch[4])]
        p_ecmwf = p_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                          (month_ecmwf == epoch[3]) | (month_ecmwf == epoch[4])]

    if epoch_text == 'summer':
        mjd_ecmwf = mjd_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2])]
        h_ecmwf = h_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2])]
        n_ecmwf = n_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2])]
        p_ecmwf = p_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2])]

    if epoch_text == 'intermediate':
        mjd_ecmwf = mjd_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                              (month_ecmwf == epoch[3])]
        h_ecmwf = h_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                          (month_ecmwf == epoch[3])]
        n_ecmwf = n_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                          (month_ecmwf == epoch[3])]
        p_ecmwf = p_ecmwf[(month_ecmwf == epoch[0]) | (month_ecmwf == epoch[1]) | (month_ecmwf == epoch[2]) |
                          (month_ecmwf == epoch[3])]

    if epoch_text == 'all':
        mjd_ecmwf = mjd_ecmwf
        h_ecmwf = h_ecmwf
        n_ecmwf = n_ecmwf
        p_ecmwf = p_ecmwf



    return mjd_ecmwf, year_ecmwf, month_ecmwf, day_ecmwf, hour_ecmwf, h_ecmwf, p_ecmwf, n_ecmwf