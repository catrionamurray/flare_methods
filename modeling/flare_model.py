from imports import *

def fit_davenport_model(time, flux, istart, istop, ipeak, buffer, small_flare, debug):
    # ADAPTED FROM FLATWRN
    # We also save the stdev for calculating start/end times
    # stdev = np.std((f - model)[regressor.inlier_mask_])

    # midflare = (time[istart] + time[istop]) / 2.
    try:
        if small_flare:
            window_mask = (time > time[max(istart - 5, 0)]) \
                          * (time < time[min(istop, len(time) - 1)] + buffer)
        else:
            window_mask = (time > time[max(istart - 10, 0)]) \
                          * (time < time[min(istop, len(time) - 1)] + buffer)
        t = time[window_mask]
        f = flux[window_mask]
        # print(max(istart - 5,0),time[max(istart - 5,0)],t[0],time[ipeak],t[ipeak - max(istart - 5,0)-1])
        # ipeak = ipeak - max(istart - 5,0)
        if debug:
            plt.plot(t, f, 'k.')
            plt.plot(time[istart:istop], flux[istart:istop], 'r.')
            plt.plot(time[ipeak], flux[ipeak], 'rx')
            plt.show(block=True)
            plt.close()

        f = f - 1  # np.nanmedian(f[0:ipeak - max(istart - 5,0)-2])
    except Exception as e:
        print(e)

    try:
        global fwhm  # I know, there is a special place in hell for this...
        if fwhm == 0:  # not defined as command-line argument
            fwhm = 1. / 24  # Selected by educated random guessing
    except NameError:
        fwhm = 1. / 24  # If calling just this function this might be handy

    # imax = np.nanargmax(f)
    # print(ipeak,np.nanargmax(flux),np.nanargmax(flux[0:ipeak+2]))
    # imax = ipeak - max(istart - 5,0)-1
    tpeak = time[ipeak]  # t[imax]
    ampl = flux[ipeak] - 1  # f[imax]#np.nanmax(f)
    if not np.isfinite(ampl):
        ampl = flux[np.nanargmax(flux)]
    fwhm = 0.25 * (t[-1] - t[0])

    pguess = (tpeak, fwhm, ampl)

    try:
        if small_flare:
            minflaret = t[0]  # t[imax-2]
            maxflaret = time[min(ipeak + 3, len(time) - 1)]  # t[imax+2]
        else:
            minflaret = t[0]
            maxflaret = time[min(ipeak + 20, len(time) - 1)]  # t[imax+20]

        popt1, pcov = curve_fit(aflare1, t, f, p0=pguess,
                                bounds=([minflaret, 0, ampl], [maxflaret, 3 * fwhm, 3 * ampl]))
    except ValueError:
        # tried to fit bad data, so just fill in with NaN's
        # shouldn't happen often
        popt1 = np.array([np.nan, np.nan, np.nan])
    except RuntimeError:
        # could not converge on a fit with aflare
        # fill with bad flag values
        popt1 = np.array([-99., -99., -99.])

    print("[PEAK FLARE T, FWHM, AMPLITUDE]")
    print("Initial guess:", pguess)
    print("Fitted result:", popt1)

    flare_t = np.linspace(np.min(t), np.max(t), int((t[-1] - t[0]) * 100000))
    flare_f = aflare1(flare_t, popt1[0], popt1[1], popt1[2])
    if debug:
        plt.plot(flare_t, flare_f, 'c')
        plt.plot(t, f, 'k.')
        plt.show()
        plt.close()

    return popt1, flare_t, flare_f