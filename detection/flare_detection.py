from utils import *
from modeling.flare_model import fit_davenport_model


def get_flare_peak(flare_start, flare_stop, peaks_old, f, t):
    peaks_new = []
    pm = 2

    for p, start, stop in zip(peaks_old, flare_start, flare_stop):
        # flareregion = f[start:stop]
        # flareregion_t = t[start:stop]
        peakregion = f[max(p - pm, start):min(p + pm, stop) + 1]
        peak = p - start  # min(p + pm, stop) - max(p - pm, start)
        st1 = max(p - pm, start) - start
        # en1 = min(p + pm, stop) - start
        # plt.plot(flareregion_t, flareregion, 'k.')
        # plt.plot(t[p], f[p], 'rx')
        # plt.show()
        # plt.close()
        try:
            while np.argmax(peakregion) != pm:
                print("Not identified the correct peak!")
                peak = np.argmax(peakregion) + st1
                peakregion = f[max(peak - pm + start, start):min(peak + pm + start, stop) + 1]
                st1 = st1 + (np.argmax(peakregion) - pm)
        except Exception as e:
            print(e)

        peak = peak + start

        # if we get the first or last point then go back to original peak
        if peak == start or peak == stop:
            peak = p

        peaks_new.append(peak)

        # plt.plot(flareregion_t, flareregion, 'k.')
        # plt.plot(t[peak], f[peak], 'rx')
        # plt.show()
        # plt.close()

    return peaks_new


def flare_start_stop(flaremask, i_peak):
    flares_start, flares_stop, flares_peak = [], [], []
    istart, istop = 0, 0
    fprev = flaremask[0]
    num_flares = len(i_peak)

    for f in range(len(flaremask)):
        if flaremask[f] == 1:
            if fprev == 0:
                # start of new flare
                istart = f

            # else:
            # continuation of flare region

            # if there's more than one flare in a flare region then split them
            # if (len(flares_start))!=num_flares:
            #     if f > (i_peak[len(flares_start)]-5):
            #         istop=f
            #         flares_start.append(istart)
            #         flares_stop.append(istop)
            #         istart = f + 1


        else:
            # not in any flare region
            if fprev == 1:
                # if end of flare
                istop = f

                # SPLIT FLARES IF MORE THAN ONE IN REGION:
                count_flares = 0
                # peaks must be sorted
                for i in range(len(i_peak)):
                    if i_peak[i] > istart and i_peak[i] < istop:
                        if count_flares > 0:
                            # if more than 1 peak is in the flare region then add new flare
                            flares_start.append(istart)
                            flares_stop.append(i_peak[i] - 2)
                            istart = i_peak[i] - 2
                        count_flares = count_flares + 1

                flares_start.append(istart)
                flares_stop.append(istop)

        fprev = flaremask[f]

    if len(flares_start) != len(flares_stop):
        flares_stop.append(f)

    return flares_start, flares_stop, i_peak


def get_flare_region(flaremask, flarepeaks):
    flares_start, flares_stop, flares_peak = [], [], []
    istart, istop, ipeak = 0, 0, 0
    fprev = flaremask[0]

    for f in range(len(flaremask)):
        if flaremask[f] == 1:
            if fprev == 0:
                istart = f
        else:
            if fprev == 1:
                istop = f
                flares_start.append(istart)
                flares_stop.append(istop)

                count_flares = 0
                for i in range(len(flarepeaks)):
                    if flarepeaks[i] > istart and flarepeaks[i] < istop:
                        flares_peak.append(flarepeaks[i])
                        if count_flares > 0:
                            flares_start.append(istart)
                            flares_stop.append(istop)
                        count_flares = count_flares + 1
        fprev = flaremask[f]

    if len(flares_start) != len(flares_stop):
        flares_stop.append(f)

    # remove flares occurring within a few dataframes of each other
    flares_start = [f for s, f in sorted(zip(flares_peak, flares_start))]
    flares_stop = [f for s, f in sorted(zip(flares_peak, flares_stop))]
    flares_peak = sorted(flares_peak)
    diff = [t - s for s, t in zip(flares_peak, flares_peak[1:])]
    todel = []
    for d in diff:
        if d <= 2:
            print("Merge Flares")
            todel.append(d)
    flares_start = np.delete(flares_start, todel)
    flares_stop = np.delete(flares_stop, todel)
    flares_peak = np.delete(flares_peak, todel)

    return flares_start, flares_stop, flares_peak


def isolate_flare_day(t, f, istart, istop, ipeak, fmask, smask, num_points_to_check_for_peak=20):
    ts, fs, split = split_multilc(t, f)
    prevsplit = 0

    if len(split)==0:
        onday_f = f
        onday_t = t
        exp = np.nanmedian([q - r for r, q in zip(onday_t, onday_t[1:])])
        onday_ind = [0,len(t)]
        onday_istart = istart
        onday_istop = istop

    for sp in range(len(split)):
        if (istart >= prevsplit) & (istart < split[sp]):
            onday_f = fs[sp]
            onday_t = ts[sp]
            exp = np.nanmedian([q - r for r, q in zip(onday_t, onday_t[1:])])
            onday_ind = [prevsplit,split[sp]]
            onday_istart = istart - prevsplit
            onday_istop = istop - prevsplit
            break
        elif sp == len(split)-1:
            onday_f = fs[-1]
            onday_t = ts[-1]
            exp = np.nanmedian([q - r for r, q in zip(onday_t, onday_t[1:])])
            onday_ind = [split[sp],len(t)]
            onday_istart = istart - split[sp]
            onday_istop = istop - split[sp]
        else:
            prevsplit = split[sp]

    onday_fmask = fmask[onday_ind[0]:onday_ind[1]]
    onday_smask = smask[onday_ind[0]:onday_ind[1]]

    # isolate flare
    flare_f = onday_f[onday_istart:onday_istop]
    flare_t = onday_t[onday_istart:onday_istop]

    peak_index = np.argmax(flare_f[:num_points_to_check_for_peak])

    return onday_t, onday_f, flare_t, flare_f, peak_index, onday_istart, onday_istop, onday_fmask, onday_smask, exp



# def qcheck_flares(t, f, istart, istop, ipeak, exp, fmask, localsize=124, globalsize=300, npoints=2, nsigma=3.0,
#                   smooth=0.1, fittype="spline"):
#     ok = True
#     do_plot = False
#     small_flare = False
#
#     # isolate night and flare region
#     t, f, flare_t, flare_f, peak_index, istart, istop, ipeak, flaremask = isolate_flare_day(t, f, istart, istop, ipeak,
#                                                                                             fmask)
#
#     masked_f = np.ma.masked_where(flaremask, f)
#     std_run_local = calculate_running_clipped_rms(t, masked_f.filled(np.nan), localsize * exp)
#     std_run_global = calculate_running_clipped_rms(t, masked_f.filled(np.nan), globalsize * exp)
#
#     # try median filter
#     medfilt_f_local = calculate_running_median(t, masked_f.filled(np.nan), localsize * exp)
#     medfilt_f_global = calculate_running_median(t, masked_f.filled(np.nan), globalsize * exp)
#
#     # ******* GET MEDIAN AND STD DEV MODELS *********
#     if fittype == "medfilt":
#         print("Fitting the lightcurve with a median filter")
#         median_model = medfilt_f_local
#         median_model_global = medfilt_f_global
#         std = std_run_local
#         std_global = std_run_global
#     elif fittype == "interp":
#         print("Fitting the lightcurve with an interpolated median filter")
#         mask = masked_f.mask
#         # ensure the last 2 data points are not masked to avoid issues with interpolation
#         mask[-2:] = False
#         mask[:2] = False
#         clean_t, clean_f, clean_medfilt, clean_medfilt_global, clean_std, clean_std_global = t[~mask], masked_f[~mask], \
#                                                                                              medfilt_f_local[~mask], \
#                                                                                              medfilt_f_global[~mask], \
#                                                                                              std_run_local[~mask], \
#                                                                                              std_run_global[~mask]
#         try:
#             # if the last/first element of the median filter is masked then set a standard fill value of 1 for interpolation
#             if clean_medfilt.mask[-1] == False and clean_medfilt.mask[0] == False:
#                 fillval = clean_medfilt[clean_medfilt.mask == False][0]
#             else:
#                 fillval = 1.0
#
#             clean_medfilt = clean_medfilt.filled(fill_value=np.ma.median(clean_medfilt))
#             clean_std = clean_std.filled(fill_value=np.ma.median(clean_std))
#
#             interpflare = interpolate.interp1d(clean_t, clean_medfilt, fill_value=fillval)
#             interp_fillgap = interpflare(t)
#             interpflare_std = interpolate.interp1d(clean_t, clean_std, fill_value=fillval)
#             interp_fillgap_std = interpflare_std(t)
#             interpflare_global = interpolate.interp1d(clean_t, clean_medfilt_global, fill_value=fillval)
#             interp_fillgap_global = interpflare_global(t)
#             interpflare_global_std = interpolate.interp1d(clean_t, clean_std_global, fill_value=fillval)
#             interp_fillgap_global_std = interpflare_global_std(t)
#             median_model = interp_fillgap
#             median_model_global = interp_fillgap_global
#             std = interp_fillgap_std
#             std_global = interp_fillgap_global_std
#
#             # TEST
#             std_run_local2 = calculate_running_clipped_rms(t, masked_f.filled(np.nan) / interp_fillgap, localsize * exp)
#             interpflare_std2 = interpolate.interp1d(clean_t, std_run_local2[~mask], fill_value=fillval)
#             interp_fillgap_std2 = interpflare_std2(t)
#         except Exception as e:
#             print(e)
#             median_model = medfilt_f_local
#             median_model_global = medfilt_f_global
#             std = std_run_local
#             std_global = std_run_global
#
#     else:
#         print("WARNING: No median fitting chosen - average median of points before and after flare region")
#         median_model = (np.nanmedian(masked_f[max(istart - 5, 0):istart + 2]) + np.nanmedian(
#             masked_f[istop - 2:min(istop + 5, len(t))])) / 2
#         median_model_global = (np.nanmedian(masked_f[max(istart - 5, 0):istart]) + np.nanmedian(
#             masked_f[min(istop + 5, len(t)):min(istop + 10, len(t))])) / 2
#
#     # *********************************
#
#     print(istop, istart, peak_index, ipeak, istop - (istart + peak_index))
#
#     # CONDITION 1: more than 2 points in flare region after flare peak:
#     if istop - (istart + peak_index) <= 2:
#         print("Flare failed Condition 1: only 1 point is in the flare region after the peak - likely a cosmic")
#         ok = False
#         flaremask[istart:istop] = 0
#
#     # CONDITION 2: more than 1 point nsigma (local) above the median:
#     above_local_thresh = flare_f - (median_model[istart:istop] + (nsigma * std[istart:istop]))
#     if len(above_local_thresh[above_local_thresh >= 0]) < npoints:
#         print(
#             "Flare failed Condition 2: only 1 point is >" + str(nsigma) + " * sigma above the median - likely a cosmic")
#         ok = False
#         flaremask[istart:istop] = 0
#
#     # CONDITION 3: more than 1 point 7*sigma (global) above the median:
#     above_global_thresh = flare_f - (median_model_global[istart:istop] + (7 * std_global[istart:istop]))
#     if len(above_global_thresh[above_global_thresh >= 0]) < npoints:
#         print("Small Flare!")
#         small_flare = True
#         boxsize = localsize
#         # if fittype == "medfilt":
#         #     print("Fitting the lightcurve with a median filter")
#         #     median_model = medfilt_f_local
#         # elif fittype == "interp":
#         #     print("Fitting the lightcurve with an interpolated median filter")
#         #     median_model = interp_fillgap
#         # else:
#         #     print("WARNING: No median fitting chosen - average median of points before and after flare region")
#         #     median_model = (np.nanmedian(masked_f[max(istart - 5, 0):istart + 2]) + np.nanmedian(
#         #         masked_f[istop - 2:min(istop + 5, len(t))])) / 2
#     else:
#         # BIG FLARE
#         print("Big Flare!")
#         boxsize = globalsize
#         median_model = median_model_global
#         std = std_global
#         # if fittype == "medfilt":
#         #     print("Fitting the lightcurve with a median filter")
#         #     median_model = medfilt_f_global
#         # elif fittype == "interp":
#         #
#         # else:
#         #     print("WARNING: No median fitting chosen - average median of points before and after flare region")
#         #     median_model = (np.nanmedian(masked_f[max(istart - 5, 0):istart]) + np.nanmedian(
#         #         masked_f[min(istop + 5 , len(t)):min(istop + 10, len(t))])) / 2
#
#     plot_area = [max(istart - 20, 0), min(istop + 20, len(t) - 1)]
#
#     x1 = std  # [plot_area[0]:plot_area[1]]
#     x2 = std_run_local  # calculate_running_rms(t, masked_f, localsize * exp)#[plot_area[0]:plot_area[1]]
#     x3 = std_run_global  # calculate_running_rms(t, masked_f, globalsize * exp)
#     # x2_2 = std_run_local2
#
#     m1 = median_model  # [plot_area[0]:plot_area[1]]
#     m2 = medfilt_f_local  # [plot_area[0]:plot_area[1]]
#     m3 = medfilt_f_global
#
#     if do_plot:
#         fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 16))
#         ax1.plot(t, f, 'k.')
#         ax1.plot(t[istart - 5:istart + 2], f[istart - 5:istart + 2], 'b.')
#
#         ax2.plot(t[plot_area[0]:plot_area[1]], f[plot_area[0]:plot_area[1]], 'k.')
#         # ax2.plot(t[plot_area[0]:plot_area[1]], m2[plot_area[0]:plot_area[1]], linestyle='--', color="lime",
#         #          label="Median Filter (local=" + str(localsize) + " points)")
#         # ax2.plot(t[plot_area[0]:plot_area[1]], m3[plot_area[0]:plot_area[1]], linestyle='--', color="darkviolet",
#         #          label="Median Filter (global=" + str(globalsize) + " points)")
#
#         if ok:
#             ax1.plot(flare_t, flare_f, 'g.', label="Vetted Flare")
#         else:
#             ax1.plot(flare_t, flare_f, 'r.', label="Unlikely Flare")
#         ax1.plot(t[ipeak], f[ipeak], marker='x', color='orange')
#
#         # ax1.plot(t, m2, linestyle='--', color="lime", label="Median Filter (local=" + str(localsize) + " points)")
#         # ax1.plot(t, m3, linestyle='--', color="darkviolet",
#         #          label="Median Filter (global=" + str(globalsize) + " points)")
#         # ax1.plot(t, m2 + (nsigma * x2), linestyle='--', color="maroon",
#         #          label=str(nsigma) + " * Running StdDev (boxsize=" + str(localsize) + " points)")
#         # ax1.plot(t, m3 + (nsigma * x3), 'r--',
#         #          label=str(nsigma) + " * Running StdDev (boxsize=" + str(globalsize) + " points)")
#         # ax1.plot(t, m3 + (7 * x3), 'r--', label="7 * Running StdDev (boxsize=" + str(globalsize) + " points)")
#
#         if fittype == "interp":
#             ax1.plot(t, interp_fillgap, linestyle='--', color="b", label="interpolate")
#             ax1.plot(t, interp_fillgap + (nsigma * interp_fillgap_std), linestyle='--', color="orange",
#                      label="interpolate + " + str(nsigma) + " * stddev")
#             ax1.plot(t, interp_fillgap + (nsigma * interp_fillgap_std2), linestyle='--', color="cyan",
#                      label="interpolate + " + str(nsigma) + " * stddev2")
#
#         ax1.plot(t, m1 + (nsigma * x1), linestyle='--', color="k",
#                  label="FINAL THRESHOLD")
#         ax1.plot(t, m1 + (3 * x1), linestyle='--', color="k",
#                  label="FINAL THRESHOLD 2")
#
#         ax1.set_ylim(bottom=np.nanmin(f) - 0.05)
#         ax1.legend(loc="lower right", fontsize='x-small', bbox_to_anchor=(1.2, 0))
#         plt.show(block=True)
#         plt.close()
#
#     # if fittype=="medfilt":
#     # plot the flare, with stellar variability/systematics corrected using a median filter
#     # f_new = f/median_model #m2
#     masked_f_new = masked_f / median_model  # m2
#     # median_new = medfilt(masked_f_new, size=boxsize)
#     # median_new = calculate_running_median(t, masked_f_new.filled(np.nan), boxsize * exp)
#     std_new = calculate_running_clipped_rms(t, masked_f_new.filled(np.nan), boxsize * exp)
#
#     # snr = ((f[istart + peak_index]-1)/np.nanmedian(std_new))
#     # print(snr)
#
#     if ok:
#         buffer = 0.5 / 24.  # days
#
#         print(istop, len(t), min(istop, len(t) - 1))
#         window_mask = (t > t[max(istart - 5, 0)]) \
#                       * (t < t[min(istop - 1, len(t) - 1)] + buffer)
#         tflare = t[window_mask]
#         fflare = f[window_mask]
#         medmodelflare = median_model[window_mask]
#
#         # model the flare on the corrected flare
#         popt, tmodel, fmodel = fit_davenport_model(t, f / median_model, istart, istop - 1, ipeak, buffer, small_flare,
#                                                    False)
#
#         if np.isnan(popt[0]):
#             print("Modelling failed")
#
#         try:
#             imod = find_nearest(t, tmodel[0])
#             imod2 = find_nearest(t, tmodel[-1])
#             # remove the median correction from the model to see the flare on the original LC
#             # if fittype=="medfilt":
#             med_interp = interpolate.interp1d(t[imod:imod2 + 1], median_model[imod:imod2 + 1],
#                                               fill_value=clean_medfilt[-1])
#             fmodel_corrected = fmodel + med_interp(tmodel)
#
#             if do_plot:
#                 fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 16))
#                 ax1.plot(t, f, 'k.')
#                 ax1.plot(tmodel, fmodel_corrected, 'b')
#                 ax2.plot(t[imod:imod2], f[imod:imod2], 'k.')
#                 ax2.plot(tmodel, fmodel_corrected, 'b')
#                 ax2.plot([t[istart + peak_index], t[istart + peak_index] + popt[1]],
#                          [(f[istart + peak_index] / 2.) + 0.5, (f[istart + peak_index] / 2.) + 0.5], 'k')
#                 ax2.axvline(t[istart + peak_index] + 4 * (popt[1] / np.log(2)), label="4 * Decay Time")
#                 ax2.legend()
#                 plt.show(block=True)
#                 plt.close()
#         #     # if fittype == "spline":
#         #     #         fmodel_corrected = fmodel + interpolate.splev(tmodel, tck, der=0)
#         except Exception as e:
#             print(e)
#     else:
#         tmodel, fmodel_corrected, fmodel, popt, medmodelflare, tflare, fflare = [], [], [], [], [], [], []
#
#     return ok, tmodel, fmodel_corrected, fmodel, popt, medmodelflare, tflare, fflare


def qcheck_flares(t,f,istart,istop,ipeak, flaremask, saturated, localsize=124,globalsize=300, npoints=2, nsigma=3.0,smooth=0.1,fittype="spline",svplt='',show=False,do_plot=True):
    ok = True
    small_flare = False
    all_t = t.copy()
    all_f = f.copy()

    all_mask = np.ma.masked_where(flaremask, f).mask
    sat_mask = np.ma.masked_where(saturated, f).mask

    # isolate night and flare region
    t, f, flare_t, flare_f, peak_index, istart, istop, flaremask, saturated, exp = isolate_flare_day(t,
                                                                                                     f,
                                                                                                     istart,
                                                                                                     istop,
                                                                                                     ipeak,
                                                                                                     flaremask,
                                                                                                     saturated)


    masked_f = np.ma.masked_where(flaremask, f)
    # medflux = (np.nanmedian(masked_f[max(istart - 5,0):istart+2]) + np.nanmedian(masked_f[istop-2:min(istop+5,len(t))]))/2
    # stdflux = np.nanstd(masked_f[max(istart - 5,0):istart+2])
    # std_run = calculate_running_rms(t,masked_f,0.05)
    std_run = calculate_running_clipped_rms(t, masked_f.filled(np.nan), 0.05)
    # std_run_local = calculate_running_rms(t, masked_f, localsize * exp)
    std_run_local = calculate_running_clipped_rms(t, masked_f.filled(np.nan), localsize * exp)
    # std_run_global = calculate_running_rms(t, masked_f, globalsize * exp)
    std_run_global = calculate_running_clipped_rms(t, masked_f.filled(np.nan), globalsize * exp)

    # try median filter
    # medfilt_f = medfilt(masked_f,size=50)
    medfilt_f = calculate_running_median(t,masked_f.filled(np.nan),50*exp)
    # medfilt_f_local = medfilt(masked_f, size=localsize)
    medfilt_f_local = calculate_running_median(t,masked_f.filled(np.nan),localsize*exp)
    # medfilt_f_global = medfilt(masked_f, size=globalsize)
    medfilt_f_global = calculate_running_median(t, masked_f.filled(np.nan), globalsize * exp)
    # medfilt_f_2 = medfilt(masked_f, size=240)
    medfilt_f_2 = calculate_running_median(t, masked_f.filled(np.nan), 240 * exp)

    # try spline fit
    mask = masked_f.mask
    clean_t, clean_f = t[~mask], masked_f[~mask]
    tck = interpolate.splrep(clean_t, clean_f,s=smooth)
    tck_global = interpolate.splrep(all_t[~all_mask], all_f[~all_mask], s=localsize)
    m = len(t)
    # s = m - math.sqrt(2 * m)
    # print(s)
    print(istop,istart,peak_index,istop - (istart + peak_index))

    # CONDITION 1: more than 2 points in flare region after flare peak:
    if istop - (istart + peak_index) <= 2:
        print("Flare failed Condition 1: only 1 point is in the flare region after the peak - likely a cosmic")
        ok = False
        flaremask[istart:istop] = 0

    # CONDITION 2: more than 1 point nsigma (local) above the median:
    above_local_thresh = flare_f - (medfilt_f_local[istart:istop] + (nsigma * std_run_local[istart:istop]))
    if len(above_local_thresh[above_local_thresh >= 0]) < npoints:
        print("Flare failed Condition 2: only 1 point is >" + str(nsigma) + " * sigma above the median - likely a cosmic")
        ok = False
        flaremask[istart:istop] = 0

    # CONDITION 3: more than 1 point 7*sigma (global) above the median:
    above_global_thresh = flare_f - (medfilt_f_global[istart:istop] + (7 * std_run_global[istart:istop]))
    if len(above_global_thresh[above_global_thresh >= 0]) < npoints:
        print("Small Flare!")
        small_flare = True
        boxsize = localsize
        if fittype == "medfilt":
            print("Fitting the lightcurve with a median filter")
            median_model = medfilt_f_local
        elif fittype == "spline":
            print("Fitting the lightcurve with a spline")
            median_model = interpolate.splev(t, tck, der=0)
        else:
            print("WARNING: No median fitting chosen - average median of points before and after flare region")
            median_model = (np.nanmedian(masked_f[max(istart - 5, 0):istart + 2]) + np.nanmedian(
                masked_f[istop - 2:min(istop + 5, len(t))])) / 2
    else:
        # BIG FLARE
        boxsize = globalsize
        if fittype == "medfilt":
            print("Fitting the lightcurve with a median filter")
            median_model = medfilt_f_global
        elif fittype == "spline":
            print("Fitting the lightcurve with a spline")
            median_model = interpolate.splev(t, tck, der=0)
        else:
            print("WARNING: No median fitting chosen - average median of points before and after flare region")
            median_model = (np.nanmedian(masked_f[max(istart - 5, 0):istart]) + np.nanmedian(
                masked_f[min(istop + 5 , len(t)):min(istop + 10, len(t))])) / 2

    # masked_f = np.ma.masked_where(flaremask,f)
    # medfilt_f = medfilt(masked_f, size=30)  # int(len_night/2.))

    plot_area = [max(istart-20,0), min(istop+20,len(t))]

    x = std_run#[plot_area[0]:plot_area[1]]
    x2 = std_run_local #calculate_running_rms(t, masked_f, localsize * exp)#[plot_area[0]:plot_area[1]]
    x3 = std_run_global #calculate_running_rms(t, masked_f, globalsize * exp)

    m1 = medfilt_f#[plot_area[0]:plot_area[1]]
    m2 = medfilt_f_local#[plot_area[0]:plot_area[1]]
    m3 = medfilt_f_global
    m4 = medfilt_f_2

    # m_interp = interp_f(t_interp)

    # t_plot = t[plot_area[0]:plot_area[1]]
    # f_plot = f[plot_area[0]:plot_area[1]]
    # masked_f_plot = masked_f[plot_area[0]:plot_area[1]]

    if do_plot:
        fig,(ax1,ax2) = plt.subplots(2,1,figsize=(16,16))
        ax1.plot(t, f, 'k.')
        ax1.plot(t[istart - 5:istart + 2],f[istart - 5:istart + 2],'b.')

        ax2.plot(t[plot_area[0]:plot_area[1]], f[plot_area[0]:plot_area[1]], 'k.')
        ax2.plot(t[plot_area[0]:plot_area[1]][saturated[plot_area[0]:plot_area[1]]==1], f[plot_area[0]:plot_area[1]][saturated[plot_area[0]:plot_area[1]]==1], 'b.')
        ax2.plot(t[plot_area[0]:plot_area[1]],m1[plot_area[0]:plot_area[1]],'c--',label="Median Filter (size=50 points)")
        ax2.plot(t[plot_area[0]:plot_area[1]],m2[plot_area[0]:plot_area[1]], linestyle='--',color="lime",label="Median Filter (local="+str(localsize)+" points)")
        ax2.plot(t[plot_area[0]:plot_area[1]],m3[plot_area[0]:plot_area[1]], linestyle='--', color="darkviolet",label="Median Filter (global="+str(globalsize)+" points)")
        ax2.plot(t[plot_area[0]:plot_area[1]],m4[plot_area[0]:plot_area[1]], linestyle='--', color="orange", label="Median Filter (size=240 points)")

        if ok:
            ax1.plot(flare_t, flare_f, 'g.', label="Vetted Flare")
        else:
            ax1.plot(flare_t, flare_f, 'r.',label="Unlikely Flare")

        ax1.plot(t[saturated==1],f[saturated==1],'c.',label="Saturated")
        ax1.plot(t,m1,'c--',label="Median Filter (size=50 points)")
        ax1.plot(t, m2, linestyle='--',color="lime",label="Median Filter (local="+str(localsize)+" points)")
        ax1.plot(t, m3, linestyle='--', color="darkviolet",label="Median Filter (global="+str(globalsize)+" points)")
        ax1.plot(t, m4, linestyle='--', color="orange", label="Median Filter (size=240 points)")

        if fittype=="spline":
            t_interp = np.linspace(min(t), max(t), num=1000, endpoint=True)
            m_interp = interpolate.splev(t_interp, tck, der=0)
            ax1.plot(t_interp, m_interp, linestyle='--', color="orange",label="Spline Fit (smooth=" + str(smooth) + ")")

            m_interp_global = interpolate.splev(t_interp, tck_global, der=0)
            ax1.plot(t_interp, m_interp_global, linestyle='--', color="orange", label="Spline Fit (smooth=" + str(localsize) + ")")

        ax1.plot(t, m1 + (nsigma * x), 'y--',label= str(nsigma) + " * Running StdDev (boxsize="+str(0.05)+"d)")
        ax1.plot(t, m2 + (nsigma * x2), linestyle='--',color="maroon",label=str(nsigma) + " * Running StdDev (boxsize="+str(localsize)+" points)")
        ax1.plot(t, m3 + (nsigma * x3), 'r--',label=str(nsigma) + " * Running StdDev (boxsize="+str(globalsize)+" points)")
        ax1.plot(t, m3 + (7 * x3), 'r--',label="7 * Running StdDev (boxsize="+str(globalsize)+" points)")
        ax1.set_ylim(bottom = np.nanmin(f)-0.05)
        ax1.legend(loc="lower right",fontsize='x-small',bbox_to_anchor=(1.2, 0))

        if show:
            plt.show()
        if svplt != "":
            plt.savefig(svplt + "_flarecheck")

        plt.close()

    # print(f[istart + peak_index])

    if fittype=="medfilt":
        # plot the flare, with stellar variability/systematics corrected using a median filter
        f_new = f/median_model #m2
        masked_f_new = masked_f/median_model #m2
        # median_new = medfilt(masked_f_new, size=boxsize)
        median_new = calculate_running_median(t, masked_f_new.filled(np.nan), boxsize * exp)
        std_new = calculate_running_clipped_rms(t, masked_f_new.filled(np.nan), boxsize * exp)

        snr = ((f[istart + peak_index]-1)/np.nanmedian(std_new))
        print(snr)

        if do_plot:
            plt.plot(t,f_new,'k.')
            plt.plot(t, median_new,'k--')
            plt.plot(t, median_new + (nsigma * std_new), linestyle='--', color="salmon")

            plt.title("Median Filter (size = " + str(boxsize) + " points)")
            if show:
                plt.show()
            if svplt != "":
                plt.savefig(svplt + "_flarecheck_medfilt")
            plt.close()


    if fittype=="spline":
        print("Number of points in LC: ",len(t))
        # plot the flare, with stellar variability/systematics corrected using spline fit
        f_new = f/interpolate.splev(t, tck, der=0)
        masked_f_new = masked_f/interpolate.splev(t, tck, der=0)
        # median_new = medfilt(masked_f_new, size=boxsize)
        median_new = calculate_running_median(t, masked_f_new.filled(np.nan), boxsize * exp)
        m_new = median_new#[plot_area[0]:plot_area[1]]
        # std_new = calculate_running_rms(t, masked_f_new, localsize * exp)#[plot_area[0]:plot_area[1]]
        std_new = calculate_running_clipped_rms(t, masked_f_new.filled(np.nan), boxsize * exp)

        snr = ((f[istart + peak_index]-1)/np.nanmedian(std_new))
        print("SNR: ",snr)

        if do_plot:
            plt.plot(t, f_new,'k.')
            plt.plot(t, m_new,'k--')
            plt.plot(t, m_new + (nsigma * std_new), linestyle='--', color="salmon")
            plt.title("Spline")

            if show:
                plt.show()
            if svplt != "":
                plt.savefig(svplt + "_flarecheck_spline")
            plt.close()

    if ok:
        buffer = 0.5 / 24.  # days

        # model the flare on the corrected flare
        popt, tmodel, fmodel = fit_davenport_model(t, f/median_model, istart, istop-1, ipeak, buffer, small_flare, debug=False)

        try:
            imod = find_nearest(t, tmodel[0])
            imod2 = find_nearest(t, tmodel[-1])
            # remove the median correction from the model to see the flare on the original LC
            if fittype=="medfilt":

                medfilt_interp = interpolate.interp1d(t[imod:imod2 + 1], medfilt_f_local[imod:imod2 + 1])

                fmodel_corrected = fmodel + medfilt_interp(tmodel)

                if do_plot:
                    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(16,16))
                    ax1.plot(t, f, 'k.')
                    ax1.plot(tmodel, fmodel_corrected, 'b')
                    ax2.plot(t[imod:imod2],f[imod:imod2],'k.')
                    ax2.plot(tmodel, fmodel_corrected, 'b')
                    ax2.plot([t[istart+peak_index],t[istart+peak_index]+popt[1]],[(f[istart+peak_index]/2.) + 0.5,(f[istart+peak_index]/2.) + 0.5],'k')
                    ax2.axvline(t[istart + peak_index] + 4*(popt[1]/np.log(2)),label="4 * Decay Time")
                    ax2.legend()
                    if show:
                        plt.show()
                    if svplt != "":
                        plt.savefig(svplt + "_flaremodel_medfilt")
                    plt.close()
            if fittype == "spline":
                fmodel_corrected = fmodel + interpolate.splev(tmodel, tck, der=0)
                if do_plot:
                    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(16,16))
                    ax1.plot(t, f, 'k.')
                    ax1.plot(tmodel, fmodel_corrected, 'b')
                    ax2.plot(t[imod:imod2],f[imod:imod2],'k.')
                    ax2.plot(tmodel, fmodel_corrected, 'b')
                    if show:
                        plt.show()
                    if svplt != "":
                        plt.savefig(svplt + "_flaremodel_spline")
                    plt.close()

        except Exception as e:
            print(e)
    else:
        tmodel, fmodel_corrected, fmodel, popt = [],[],[],[]

    return ok, tmodel, fmodel_corrected, fmodel, popt