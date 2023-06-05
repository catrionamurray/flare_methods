#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 16:55:22 2020

@author: fl386
"""

from imports import *

# divide light curve into nights
def daynight(tact):
    newday = (np.where((tact - np.roll(tact, 1)) > 0.5)[0])
    nd = np.zeros((np.size(newday) + 2))
    nd[1:(np.size(newday) + 1)] = newday
    nd[(np.size(newday) + 1)] = np.size(tact)
    nd = nd.astype(int)
    return nd

# plot more conveniently with t_plot (no big gaps), ndd (position of beginning of night)
# and xticks
def getpl(time, nd, zeronr):
    t_plot = []
    ndd = []
    xticks = []
    for ie in range(np.size(nd) - 1):

        if ie == 0:
            tlast = 0
            tadd = 0
        else:
            tlast = t_plot[-1]
            tadd = time[nd[ie]] - time[nd[ie] - 1]

        if tadd > 10:
            kkjdd = 0.5
        else:
            kkjdd = 0.05
        t_plot.extend(time[nd[ie]:nd[ie + 1]] - time[nd[ie]] + tlast + kkjdd)
        xticks.extend([np.round(time[nd[ie]] - zeronr, 1)])
        ndd.extend([tlast + kkjdd])

    t_plot = np.array(t_plot)
    ndd = np.asarray(ndd)
    xticks = np.asarray(xticks)
    return t_plot, ndd, xticks


"""
#this is how I defined it in the paper. gives the same results.
def fp(yval,rmsh,threshval):
    yylo=(abs(yval - np.roll(yval,2))-abs(np.roll(yval,-1)-yval))/rmsh
    yyl=(2.*yval - np.roll(yval,2) - np.roll(yval,-3))/rmsh
    yyl[yyl<0]=0
    return np.where(yylo*yyl>threshval)[0]
"""


def decay_two(xhh, aa):
    try:
        x = aa[0] + aa[1] * np.exp(-(xhh - np.nanmin(xhh)) / aa[2]) + aa[3] * np.exp(-(xhh - np.min(xhh)) / aa[4])
    except Exception as e:
        print(e)
    return x#aa[0] + aa[1] * np.exp(-(xhh - np.nanmin(xhh)) / aa[2]) #+ aa[3] * np.exp(-(xhh - np.min(xhh)) / aa[4])

def decay(xhh, aa):
    try:
        x = aa[0] + aa[1] * np.exp(-(xhh - np.nanmin(xhh)) / aa[2]) #+ aa[3] * np.exp(-(xhh - np.min(xhh)) / aa[4])
    except Exception as e:
        print(e)
    return x#aa[0] + aa[1] * np.exp(-(xhh - np.nanmin(xhh)) / aa[2]) #+ aa[3] * np.exp(-(xhh - np.min(xhh)) / aa[4])


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def clipped_std(y):
    return np.ma.std(sigma_clip(y, sigma_upper=3, sigma_lower=8))

def running_box(x,y,boxsize,operation):

    if operation == "std":
        op = np.nanstd
    elif operation == 'median':
        op = np.nanmedian
    elif operation == 'mean':
        op = np.nanmean
    elif operation == "clipped_std":
        op = clipped_std
    else:
        print("No running function selected - choose std, median or mean.")
        return np.nan

    l = len(x)
    boxsize = boxsize / 2
    dy = np.zeros((np.size(x)))

    for i in range(0, l):
        s = x[i]
        box_s = find_nearest(x, s - boxsize)
        box_e = find_nearest(x, s + boxsize)

        # calculate [OPERATION] for the current window
        try:
            y_ext = y[box_s:box_e + 1]
            d = op(y_ext)
            if not np.isnan(d):
                dy[i] = d

        except Exception as e:
            print(e)

    dy = np.ma.masked_where(dy == 0, dy)
    # dy = np.ma.masked_where(np.isnan(dy), dy)
    dy = np.ma.masked_invalid(dy)
    # dy = dy.filled(np.nanmedian(dy))
    # print("RUNNING " + operation)
    # print("FILL VALUE = " + str(np.nanmedian(dy)))
    return dy

# def running_box(x,y,boxsize,operation):
#     if operation == 'median':
#         op = np.nanmedian
#     elif operation == "std":
#         op = clipped_std
#
#     l = len(x)
#     boxsize = boxsize / 2
#     dy = np.zeros((np.size(x)))
#
#     for i in range(0, l):
#         # s = x[i]
#         box_s = int(max(0,i-boxsize))#find_nearest(x, s - boxsize)
#         box_e = int(min(l,i+boxsize))#find_nearest(x, s + boxsize)
#         # calculate [OPERATION] for the current window
#         try:
#             y_ext = y[box_s:box_e + 1]
#             d = op(y_ext)
#             if not np.isnan(d):
#                 dy[i] = d
#
#             # if np.isnan(d):
#             #     print("confuse")
#
#         except Exception as e:
#             # print(e)
#             pass
#
#     dy = np.ma.masked_where(dy == 0, dy)
#     dy = np.ma.masked_invalid(dy)
#     # dy = dy.filled(np.nanmedian(dy))
#
#     # prev=0
#     # next=0
#     #
#     # for d in range(len(dy)):
#     #     if dy[d]==0 or np.isnan(dy[d]):
#     #         if d > 0:
#     #             prev =d-1
#     #         if d < len(dy):
#     #             next= d+1
#     #
#     #         try:
#     #             while dy[prev] == 0 or np.isnan(dy[prev]):
#     #                 prev = prev - 1
#     #         except:
#     #             prev = 0
#     #
#     #         try:
#     #             while dy[next] == 0 or np.isnan(dy[next]):
#     #                 next = next + 1
#     #         except:
#     #             next = 0
#     #
#     #         if next==0 and prev==0:
#     #             dy[d] = np.nanmedian(dy)
#     #         elif next==0:
#     #             dy[d] =dy[prev]
#     #         elif prev==0:
#     #             dy[d] = dy[next]
#     #         else:
#     #             dy[d] = (dy[prev] + dy[next]/2.)
#
#     return dy

def split_multilc(t, y):
    start = t[0]
    t = np.array(t)
    nextdays=t[np.absolute(t-start)>0.5]
    split=[]
    while nextdays!=[]:
        start = nextdays[0]
        ind_st = np.where(t==start)[0][0]
        split.append(ind_st)
        time = t[ind_st:]
        nextdays = time[np.absolute(time-start) > 0.5]

    times=np.split(t,split)
    ys=np.split(y,split)

    return times, ys, split

def isolate_flare_day(t,f,kss):
    ts, fs, split = split_multilc(t, f)
    prevsplit = 0
    # to_subtract=0

    try:

        if len(split)==0:
            onday_f = f
            onday_t = t
            to_subtract = 0
            # onday_istart = kss
            # onday_istop = kss

        for sp in range(len(split)):
            if (kss >= prevsplit) & (kss < split[sp]):
                onday_f = fs[sp]
                onday_t = ts[sp]
                # onday_istart = istart - prevsplit
                # onday_istop = istop - prevsplit
                to_subtract = prevsplit
                # onday_kss = kss-prevsplit
                break
            elif sp == len(split)-1:
                onday_f = fs[-1]
                onday_t = ts[-1]
                # onday_istart = istart - split[sp]
                # onday_istop = istop - split[sp]
                # onday_kss = kss - split[sp]
                to_subtract = split[sp]
            else:
                prevsplit = split[sp]

        # print(onday_t[0])
    except Exception as e:
        print(e)

    return onday_t,onday_f, to_subtract



def geterr20(der, ndimage, kj):
    enr = 61
    av = int(enr / 2.)
    drr = der[ndimage[kj]:(ndimage[kj + 1])]
    if (ndimage[kj + 1] - ndimage[kj] <= enr) and ndimage[kj + 1] - ndimage[kj] > 2:
        enr = (ndimage[kj + 1] - ndimage[kj])
        res1 = np.ones((np.size(drr))) * np.ma.std(sigma_clip(drr, sigma_upper=6, sigma_lower=10))

    elif (ndimage[kj + 1] - ndimage[kj] <= 2):
        res1 = np.ones((ndimage[kj + 1] - ndimage[kj])) * 1.
    else:
        midval = np.asarray([np.ma.std(sigma_clip(drr[i:i + enr], sigma_upper=6, sigma_lower=10)) for i in
                             range(np.size(drr) - 2 * av)])

        strtval = np.ones((av))
        endval = np.ones((av))
        # higher rms expected for start/end of night. -> recompute with decreasing box size.
        for j in range(av):
            strtval[av - 1 - j] = np.ma.std(sigma_clip(drr[:(enr - 1 - j)], sigma_upper=6, sigma_lower=10))
            endval[-av + j] = np.ma.std(sigma_clip(drr[-enr + 1 + j:], sigma_upper=6, sigma_lower=10))
        res1 = np.hstack((strtval, midval, endval))
    return res1


def findflares_sub(threshval, k, lc, xval, ndimage, fluxdiff, fluxdiff2, rmsh,rmsh_test, exp, run,boxsize_local,boxsize_global,debug=False):
    # removed condition  np.max(testlc)>(medianh+2.*stdh) or np.max(testlc)>1.02*medianh:

    # threshval: threshold for condition 1
    # k: number of sigma for detection?
    # lc: light curve
    # xval: jd
    # ndimage: where do new nights start?
    # rmsh box rms
    # t_axis: other time values. ordered like jd.
    # run: run cutflares 2 to 3 times. run = number of run

    # ------------------------------------------PART 1: FIND FLARES

    # # fj - fj-1
    # fluxdiff = (lc - np.roll(lc, 1))
    # # fj - fj-2
    # fluxdiff2 = (lc - np.roll(lc, 2))

    # second half of condition 1
    # yylo = (abs(fluxdiff2) - abs(np.roll(fluxdiff, -1))) / rmsh
    yylo = (abs(lc - np.roll(lc, 2)) - abs(lc - np.roll(lc, -1))) / rmsh

    # first half of condition 1
    yyl = (fluxdiff2 - np.roll(fluxdiff, -1) - np.roll(fluxdiff, -2) - np.roll(fluxdiff, -3)) / rmsh
    yyl2 = ((2*lc) - np.roll(lc, 2) - np.roll(lc, -3))/ rmsh
    yyl3 = ((2 * lc) - np.roll(lc, 2) - np.roll(lc, -2))/ rmsh
    yyl4 = ((2 * lc) - np.roll(lc, 3) - np.roll(lc, -3)) / rmsh
    # condition 2
    yyl[yyl < 0] = 0

    yyl2[yyl2 < 0] = 0
    yyl3[yyl3 < 0] = 0
    yyl4[yyl4 < 0] = 0

    # stdval = np.ma.std(sigma_clip(lc, sigma=2))
    # stdval = max(stdval, 0.002)

    # find flares
    flareposw = np.where(yylo * yyl > threshval)[0]

    if debug:
        plt.figure(figsize=(24,6))
        plt.plot(xval, yylo*yyl, 'k.')
        plt.plot(xval, yylo * yyl2, 'r.')
        plt.ylim(bottom=0,top=10)
        # plt.xlim(633.5,635)
        plt.plot(xval[flareposw], (yylo * yyl)[flareposw], 'gx')
        plt.grid(True)
        plt.show()
        plt.close()

        # plt.plot(xval, (yylo*yyl)*rmsh/rmsh_test, 'k.')
        # plt.ylim(bottom=0,top=10)
        # plt.xlim(633.5,635)
        # plt.show()
        # plt.close()

        # plt.plot(xval, yylo*yyl, 'k.')
        # plt.plot(xval, yylo * yyl3, 'g.')
        # plt.ylim(bottom=0, top=10)
        # plt.xlim(633.5, 635)
        # plt.plot(xval[flareposw], (yylo * yyl)[flareposw], 'gx')
        # plt.show()
        # plt.close()
        #
        # plt.plot(xval, yylo*yyl, 'k.')
        # plt.plot(xval, yylo * yyl4, 'b.')
        # plt.ylim(bottom=0, top=10)
        # plt.xlim(633.5, 635)
        # plt.plot(xval[flareposw], (yylo * yyl)[flareposw], 'gx')
        # plt.show()
        # plt.close()


    # ------------------------------------------PART 2: CLEAN FLARE POSITIONS. e.g. REMOVE "FLARES" IN NOISE

    # exclude some flare candidates
    flarepos = []
    flarethresh = []
    flareend,flarestart,flaremed,flaretosub,flaredayx,flaredayy = [],[],[],[],[],[]

    cnlr = 0
    lcsize = np.size(lc)
    to_subtract = 0
    prev_kss = 0
    for kss in flareposw:
        # print xval[kss]
        cnlr += 1
        plusval = 4

        if cnlr == 1 or xval[kss] - xval[prev_kss] > 0.5:
            onday_x, onday_lc, to_subtract = isolate_flare_day(xval, lc, kss)
            onday_lcsize = np.size(onday_lc)

        onday_kss = kss - to_subtract
        prev_kss = kss

        # plt.plot(xval[kss-50],lc[kss+50],'k.')
        # plt.plot(onday_x,onday_lc,'k.')
        # plt.plot(xval[kss],lc[kss],'r.',markersize=12)
        # plt.show()
        # plt.close()

        if onday_kss < 1:
            continue

        # flare at end of lc? -> not enough data points for fit. -> will be removed in cosmic detection
        if (onday_kss + plusval) > onday_lcsize - 1:
            continue

        # flare next to other flare?
        if cnlr > 1:
            if kss < (flareposw[cnlr - 2] + 3) and (xval[kss] - xval[(flareposw[cnlr - 2] + 3)]) < 0.5:
                continue
            # if kss - flarepos[-1] < 60:
            #     continue

        # flare right at the end of a night (/gap)?
        if (xval[kss + 3] - xval[kss]) > 0.2:
            continue

        # include some data before flare starts, allow at least 2 points before and after flare to calculate median/rms
        minusval = 10#3 #2

        # except if flare happens in the first few data points
        if onday_kss < minusval + 2:#3:
            # minusval = onday_kss-2
            minusval = onday_kss


        # significant?
        if onday_kss + 62 <= onday_lcsize:
            plusval = 60
        else:
            plusval = onday_lcsize - onday_kss - 2

        # extract small LC section around flare and sigmaclip
        # *** Cat - I don't see why you wouldn't take the average of the regions before and after...
        # or do a median filter, to get the average LC value. If the flare is large then clipping won't
        # work well ***
        xpoints = 5
        npoints =2

        # if (kss + plusval + xpoints <= lcsize) and (kss - minusval > xpoints):
        #     testlc1 = lc[kss - minusval-xpoints:kss - minusval]
        #     testlc2 = lc[kss + plusval:kss + plusval + xpoints]
        #     testlc = np.concatenate((testlc1,testlc2))
        #     medianh1 = np.ma.median(testlc)

        # boxsize_local = 80
        # boxsize_global = 320
        flaremask = np.zeros(onday_lcsize)
        # flaremask[onday_kss - minusval:onday_kss + 10] = 1
        flaremask[onday_kss - minusval:onday_kss + plusval] = 1
        if len(flarepos)>0:
            for f in flarepos:
                if xval[kss] - xval[f] < 0.5:
                    flaremask[f-to_subtract:f-to_subtract+10]=1

        # mask all flars in LC
        masked_lc = np.ma.masked_where(flaremask, onday_lc)

        # calculate running median - for small and large boxsize
        med_local = running_box(onday_x, masked_lc.filled(np.nan), boxsize_local*exp, operation="median")
        med_global = running_box(onday_x, masked_lc.filled(np.nan), boxsize_global*exp, operation="median")
        medianh = np.ma.median(med_local[onday_kss - minusval:onday_kss + plusval])
        medianh_global =np.nanmedian(med_global[onday_kss - minusval:onday_kss + plusval])

        # calculate running std dev - for small and large boxsize
        # onday_rmsh= running_box(onday_x, onday_lc, boxsize_local,operation="std")
        # onday_rmsh_global = running_box(onday_x, onday_lc, boxsize_global,operation="std")
        onday_rmsh= running_box(onday_x, masked_lc.filled(np.nan), boxsize_local*exp,operation="std")
        onday_rmsh_global = running_box(onday_x, masked_lc.filled(np.nan), boxsize_global*exp,operation="std")


        # TEST INTERPOLATION

        try:
            # print("Fitting the lightcurve with an interpolated median filter")
            mask = masked_lc.mask
            # ensure the last 2 data points are not masked to avoid issues with interpolation
            mask[-2:] = False
            mask[:2] = False
            clean_t, clean_f, clean_medfilt, clean_medfilt_global, clean_std, clean_std_global = onday_x[~mask], masked_lc[~mask], \
                                                                                                 med_local[~mask], \
                                                                                                 med_global[~mask], \
                                                                                                 onday_rmsh[~mask], \
                                                                                                 onday_rmsh_global[~mask]
            # if the last/first element of the median filter is masked then set a standard fill value of 1 for interpolation
            if clean_medfilt.mask[-1] == False and clean_medfilt.mask[0] == False:
                fillval = clean_medfilt[clean_medfilt.mask == False][0]
            else:
                fillval = 1.0

            clean_medfilt = clean_medfilt.filled(fill_value=np.ma.median(clean_medfilt))
            clean_std = clean_std.filled(fill_value=np.ma.median(clean_std))

            interpflare = interpolate.interp1d(clean_t, clean_medfilt, fill_value=fillval)
            interp_fillgap = interpflare(onday_x)
            interpflare_std = interpolate.interp1d(clean_t, clean_std, fill_value=fillval)
            interp_fillgap_std = interpflare_std(onday_x)
            interpflare_global = interpolate.interp1d(clean_t, clean_medfilt_global, fill_value=fillval)
            interp_fillgap_global = interpflare_global(onday_x)
            interpflare_global_std = interpolate.interp1d(clean_t, clean_std_global, fill_value=fillval)
            interp_fillgap_global_std = interpflare_global_std(onday_x)
            median_model = interp_fillgap
            median_model_global = interp_fillgap_global
            std = interp_fillgap_std
            std_global = interp_fillgap_global_std
        except Exception as e:
            print(e)
            median_model=medianh
            median_model_global = medianh_global
            std = onday_rmsh
            std_global = onday_rmsh_global

        testlc =  lc[kss - minusval:kss + plusval]
        testt = xval[kss - minusval:kss + plusval]
        mlocal = median_model[onday_kss - minusval:onday_kss + plusval]
        slocal = std[onday_kss - minusval:onday_kss + plusval]
        mglobal = median_model_global[onday_kss - minusval:onday_kss + plusval]
        sglobal = std_global[onday_kss - minusval:onday_kss + plusval]

        testlc_clipped = sigma_clip(testlc, sigma_upper=3, sigma_lower=8)
        medianh2 = np.ma.median(testlc_clipped)

        # THRESHOLD TO REMOVE COSMICS
        above_local_thresh = testlc[:minusval+10] - (mlocal[:minusval+10] + k * slocal[:minusval+10])
        # above_local_thresh = above_local_thresh[~np.isnan(above_local_thresh)]
        above_local_thresh = above_local_thresh[above_local_thresh.mask == False]
        if len(above_local_thresh[above_local_thresh >= 0]) < npoints:
            continue

        above_global_thresh =  testlc[:minusval+10] - (mglobal[:minusval+10] + 7 * sglobal[:minusval+10])
        # above_global_thresh = above_global_thresh[~np.isnan(above_global_thresh)]
        above_global_thresh = above_global_thresh[above_global_thresh.mask == False]

        if len(above_global_thresh[above_global_thresh >= 0]) >= npoints:
            # BIG FLARE
            print("BIG FLARE")
            median_final = median_model_global[onday_kss]#medianh_global
            med = median_model_global
            rms_final = std_global[onday_kss]
            rms = std_global
            # bigflare = True
            # end = kss + int(7200 / (exp * 24 * 3600))
            # if onday_kss = kss - to_subtract

            # define start and end of flare region, larger for a big flare
            if onday_kss + int(7200 / (exp * 24 * 3600)) <= onday_lcsize - 2:
                end = kss + int(7200 / (exp * 24 * 3600))
            else:
                end = kss + onday_lcsize - onday_kss - 2

            if onday_kss < 12:
                srt = to_subtract + 2
            else:
                srt = kss - 10
        else:
            # SMALL FLARE
            print("SMALL FLARE")
            median_final = median_model[onday_kss]#medianh[onday_kss]
            med = median_model
            rms_final = std[onday_kss]
            rms = std
            # bigflare = False
            # end = kss + int(5400 / (exp * 24 * 3600))
            # define start and end of flare region, smaller for a small flare
            if onday_kss + int(5400 / (exp * 24 * 3600)) <= onday_lcsize - 2:
                end = kss + int(5400 / (exp * 24 * 3600))
            else:
                end = kss + onday_lcsize - onday_kss - 2

            srt = kss-minusval

        if debug:
            # plt.plot(onday_x,onday_lc,'r.')
            # plt.plot(onday_x, masked_lc, 'k.')
            # plt.show()
            # plt.close()
        #
            plt.plot(onday_x, onday_lc, 'c.')
            plt.plot(testt, testlc, 'b.')
            # plt.plot(testt[:10], testlc[:10], 'k.')
            plt.plot(xval[kss], lc[kss], 'r.')
            # plt.plot(xval[kss - minusval:kss + plusval], medianh2 + rmsh[kss - minusval:kss + plusval], 'yellow')
            plt.plot(onday_x, med_global, color='b',label="Median Global Old")
            plt.plot(onday_x, med_local, color='green',label="Median Local Old")
            plt.plot(onday_x, median_model, color="k", label="Median Local New")
            plt.plot(onday_x, median_model_global, color="k", label="Median Global New")
            plt.plot(onday_x, med, color="magenta",label="Median")
            plt.plot(onday_x, med + (k*rms), color="r",label="Median + "+str(k)+" * RMS")
            plt.plot(onday_x[onday_kss - minusval:onday_kss + plusval],
                     med[onday_kss - minusval:onday_kss + plusval] + k * rms[onday_kss - minusval:onday_kss + plusval],
                     'cyan')
            # plt.plot(onday_x,med_global,color='b')
            plt.axhline(median_final + (k * rms_final), color="green", linestyle="--",label="Final Thresh")
            plt.axhline(median_final + ((k-1) * rms_final), color="green", linestyle="--")
            plt.axhline(medianh + (k * rmsh[kss]), color='orange')
            plt.legend()
            plt.show()
            plt.close()



        if len(testlc)<=3:
            continue

        # testlc = lc[kss - minusval:kss + plusval]
        # testlc_clipped = sigma_clip(testlc, sigma_upper=3, sigma_lower=8)
        # medianh = np.ma.median(testlc_clipped)

        # exclusion criterion: if the max of the subsection > median + k*sigma
        # if np.max(testlc[:10]) > (medianh + k * rmsh[kss]):

        if np.nanmax(testlc[:minusval+10]) > (median_final + k * rms_final):

            overlap = False
            # for i in range(len(flarepos)):
            #     if flarepos[i] > srt and flarepos[i] < end:
            #         if flarestart[i] > srt:
            #             flarestart[i] = srt
            #         if flareend[i] < end:
            #             flareend[i] = end
            #         overlap = True

            if overlap == False:
                flarepos.extend([kss])
                flarethresh.extend([(yylo * yyl)[kss]])
                # flaretypes.extend([bigflare])
                flareend.extend([end])
                flarestart.extend([srt])
                flaremed.extend([med])
                flaretosub.extend([to_subtract])
                flaredayx.extend([onday_x])
                flaredayy.extend([onday_lc])

            # if debug:
            #     # debug - check plot
            #     plt.figure(figsize=(10,6))
            #     plt.plot(onday_x, onday_lc, 'c.')
            #     plt.plot(testt, testlc, 'b.')
            #     plt.plot(testt[:10], testlc[:10], 'k.')
            #     plt.plot(xval[kss], lc[kss], 'r.')
            #     plt.plot(xval[kss - minusval:kss + plusval], medianh2 + (k*rmsh[kss - minusval:kss + plusval]), 'yellow')
            #     plt.axhline(median_final + k * rms_final, color="green", linestyle="--",label="Median + "+str(k)+" * RMS (at flare)")
            #     plt.plot(onday_x, med + k * rms, color="k", label="Median + " + str(k) + " * RMS")
            #     plt.plot(onday_x[onday_kss - minusval:onday_kss + plusval],
            #              med[onday_kss - minusval:onday_kss + plusval] + (k*rms[onday_kss - minusval:onday_kss + plusval]),
            #              'cyan')
            #     plt.plot(onday_x, med, color="k", label="Median")
            #     plt.axhline(medianh + k * rmsh[kss], color='orange')
            #     plt.legend()
            #     plt.show()
            #
            #     plt.close()
        # else:
        #     print("excluded")

    flarepos = np.asarray(flarepos)
    flarethresh = np.asarray(flarethresh)
    print(np.size(flarepos),"flare candidates")
    print(flarethresh)
    print(flarepos)
    print(flarestart,flareend)

    flare_start = []
    flare_thresh = []
    decaytime = []
    i_peak = []
    cntr = 0

    if np.size(flarepos) == 0:
        # print("TEST")
        return np.array([]), np.array([]),np.array([]),np.array([])

    # ------------------------------------------PART 3: FIT FLARES

    if np.size(flarepos) >= 1:

        for kj in flarepos:

            # minusval = 3
            # if kj < 3:
            #     minusval = kj
            # if kj-flaretosub[cntr] < 3:
            #     minusval = kj-flaretosub[cntr]

            # find data points around flare. srt until end.

            # srt = kj - minusval
            srt = flarestart[cntr]
            # end = kj + 60
            # T = 60*40 = 2400s for TRAPPIST
            # n = T/te = 2400s/te for SSO
            # n = 2400s/30 = 80 points
            # n = 2400s/20 = 120 points


            # end = kj + int(2400/(exp*24*3600))
            # if flaretypes[cntr] == False:
            #     end = kj + int(5400 / (exp * 24 * 3600))
            # else:
            #     end = kj + int(7200 / (exp * 24 * 3600))
            end = flareend[cntr]
            # print(int(2400/(exp*24*3600)))

            if end > np.size(xval):
                end = np.size(xval) - 2

            # fit across two days? (we need a certain nr of data points for the fit.
            while xval[end] - xval[srt] > 0.5 and end - srt >= 5:
                end -= 1

            overlap = False
            # if flare overlaps with another flare...
            for i in flarepos:
                if i!=kj:
                    if i > srt and i < end:
                        end = i - 2#minusval
                        overlap = True
                    # elif i > end:
                    #     break

            if srt >= end:
                print("Flare of no length!")
                cntr = cntr + 1
                continue

            medsrt = srt - flaretosub[cntr]
            medend = end - flaretosub[cntr]
            yp_med = flaremed[cntr][medsrt:medend]

            # xp: preliminary flare area
            xp = xval[srt:end]
            yp_undetrend = lc[srt:end]
            yp = np.ma.divide(yp_undetrend,yp_med)

            # find flare peak
            peakindex = (np.where(yp == np.nanmax(yp[:min(10, end - srt)]))[0])[0]
            # i_peak.append(peakindex + srt)

            if peakindex<2:
                if medsrt >=2:
                    srt=srt-2
                    medsrt = srt - flaretosub[cntr]
                    medend = end - flaretosub[cntr]
                    yp_med = flaremed[cntr][medsrt:medend]

                    xp = xval[srt:end]
                    yp_undetrend = lc[srt:end]
                    yp = np.ma.divide(yp_undetrend, yp_med)
                    peakindex = (np.where(yp == np.nanmax(yp[:min(10, end - srt)]))[0])[0]



            # data points for exponential decay fit
            x_decayzone = xp[peakindex:]
            y_decayzone = yp[peakindex:]

            # time gap before peak?
            if xp[peakindex] - xp[peakindex - 2] > 0.049 or peakindex < 2:
                flarezone = xp[peakindex:]
                ctplus = 0.
            else:
                flarezone = xp[(peakindex - 2):]
                ctplus = xp[peakindex] - xp[peakindex - 2]

            # higher weight for data points affected by flare to improve fit.
            weight = np.ones((np.size(x_decayzone)))
            weight2 = np.ones((5)) * 3.

            if np.size(x_decayzone) > 5:
                lkl = 5
            else:
                lkl = np.size(x_decayzone)
            weight[:lkl] = weight2[:lkl]

            # ignore low flux data points. don't use if variability is high. not the case in trappist.
            # weight[y_decayzone < 0.98] = 0.

            # don't fit if nr of dp < 3. will be removed in cosmic removal part.
            if np.size(y_decayzone) < 3:
                print("Decay time is too short")
                cntr = cntr + 1
                continue

            # fit exponential decay
            # resk = leastsq(lambda a: (y_decayzone-decay(x_decayzone,a))*weight, p0)
            try:
                # initial guess for parameters
                p0 = (1.2, 1., 0.001, 0.5, 0.01)

                print("Initial guess for parameters: ",p0)
                if np.count_nonzero(np.isnan(y_decayzone)) > 0 :
                    print(np.count_nonzero(np.isnan(y_decayzone)))
                    print(type(y_decayzone))
                    print(type(y_decayzone.filled(np.nan)))
                    print(np.count_nonzero(np.isnan(y_decayzone.filled(np.nan))))

                    ind_todelete = np.argwhere(np.isnan(y_decayzone.filled(np.nan)))
                    # print(ind_todelete)
                    xd = np.delete(x_decayzone,ind_todelete)
                    yd = np.delete(y_decayzone, ind_todelete)

                    # higher weight for data points affected by flare to improve fit.
                    w = np.ones((np.size(xd)))
                    w2 = np.ones((5)) * 3.

                    if np.size(xd) > 5:
                        lkl = 5
                    else:
                        lkl = np.size(xd)
                    w[:lkl] = w2[:lkl]

                else:
                    xd = x_decayzone
                    yd = y_decayzone.filled(np.nan)
                    w=weight

                resk = (least_squares(lambda a: (yd - decay_two(xd, a)) * w, p0,
                                  bounds=([np.array([0, 0, 0, 0, 0]), np.array([np.inf, np.inf, np.inf, np.inf, np.inf])]))).x
                y_model = resk[0] + \
                          resk[1] * np.exp(-(x_decayzone - np.nanmin(x_decayzone)) / resk[2]) + \
                          resk[3] * np.exp(-(x_decayzone - np.nanmin(x_decayzone)) / resk[4])

                print("Final guess for parameters: ",resk)

                # if debug:
                #     # plt.plot(xp, yp_undetrend, 'r.')
                #     plt.plot(xp, yp, 'k.')
                #     plt.plot(x_decayzone, y_model, 'b')
                #     plt.show()
                #     plt.close()

            except Exception as e:
                print(e)
                print("Error with Florian Flare code")
                cntr = cntr + 1
                continue
            # print(resk)

            # characteristic decay time
            # ctime = resk[2]#
            # ctime = np.nanmean([resk[2],resk[4]])
            ctime = np.max([resk[2], resk[4]])
            if ctime > 1:
                ctime = resk[2]

            if  xval[kj] + ctime > flaredayx[cntr][-1]:
                ctime = (flaredayx[cntr][-1]-xval[kj])/4.


            # if debug:
            #     plt.plot(xp,yp_undetrend,'r.')
            #     plt.plot(xp,yp,'k.')
            #     plt.plot(x_decayzone,y_model,'b')
            #     # plt.axvline(flarezone[0] + ctime * 4 + ctplus)
            #     if ctime < 0.1:
            #         plt.axvline(flarezone[0] + ctime * 4 + ctplus,color='red')
            #     if resk[2] < 0.1:
            #         plt.axvline(flarezone[0] + resk[2] * 4 + ctplus, color='yellow')
            #     # plt.ylim(0.98,yp[peakindex]+0.01)
            #     plt.show()
            #     plt.close()

            # ignore if decay time very long. will be corrected by GP.
            if ctime > 0.1:#0.05:
                # ctime = 0
                print("Decay time is too long")
                if debug:
                    plt.plot(xp,yp_undetrend,'r.')
                    plt.plot(xp,yp,'k.')
                    plt.plot(x_decayzone,y_model,'b')
                    # plt.axvline(flarezone[0] + ctime * 4 + ctplus)
                    if ctime < 0.1:
                        plt.axvline(flarezone[0] + ctime * 4 + ctplus,color='red')
                    if resk[2] < 0.1:
                        plt.axvline(flarezone[0] + resk[2] * 4 + ctplus, color='yellow')
                    # plt.ylim(0.98,yp[peakindex]+0.01)
                    plt.show()
                    plt.close()
                cntr = cntr + 1
                continue

            flare_start.extend([flarezone[0]])
            i_peak.append(peakindex + srt)
            flare_thresh.extend([flarethresh[cntr]])
            # does the flare end before the last data point of a night?
            if ((flarezone[0] < xval[ndimage[:-1]]) & ((flarezone[0] + ctime * 4 + ctplus) > xval[ndimage[:-1]])).any():
                dy = (np.where(
                    (xval[ndimage[:-1]] > flarezone[0]) & (xval[ndimage[:-1]] < (flarezone[0] + ctime * 4 + ctplus)))[
                    0])[0]
                decaytime.extend([xval[ndimage[dy]] - flarezone[0]])
            else:
                decaytime.extend([ctime * 4 + ctplus])

            cntr = cntr + 1
        # print(flarethresh)
        return np.asarray(flare_start), np.asarray(decaytime), np.asarray(i_peak), flarethresh


def findflares(xval, yval,thresh,sigma,debug=False):
    # find beginning of "nights" (/light curve parts which are more than half a day apart from
    # other light curve parts)
    nd = daynight(xval)
    i_peak = []
    i_thresh = []
    exp = xval[1]-xval[0]
    print("EXPOSURE: " + str(24 * 60 * 60 * exp) + " seconds")
    print("Threshold = ",thresh,", Sigma = ",sigma)

    # normalise the light curves for each night
    for ie in range(np.size(nd) - 1):
        yval[nd[ie]:nd[ie + 1]] /= np.nanmedian(yval[nd[ie]:nd[ie + 1]])

    # fj - fj-1
    fluxdiff = (yval - np.roll(yval, 1))
    # fj - fj-2
    fluxdiff2 = (yval - np.roll(yval, 2))

    geterrf = partial(geterr20, fluxdiff, nd)
    # p = Pool(1)

    # rmshp = p.map(geterrf, range(len(nd) - 1))
    rmshp = []
    for r in range(len(nd) - 1):
        rmshp.append(geterrf(r))
    rmsh = np.hstack(np.array(rmshp))

    # replaceval = np.nanmedian(yval) - 0.001
    replaceval = np.nan
    if (yval == replaceval).any():
        print("change replaceval")

    print("start run1")

    rmsh_test = running_box(xval, yval, 80 * exp, operation="std")

    rss1 = findflares_sub(thresh, sigma, yval, xval, nd, fluxdiff, fluxdiff2, rmsh,rmsh_test, exp, 1,80,320,debug)
    # print(rss1)
    if np.size(rss1[0]) != 0:
        i_peak = rss1[2]
        i_thresh = rss1[3]

    print("done run1")

    # if flarefound==1: normal lc. if flarefound==0: affected by flare.
    flarefound = np.zeros((np.size(xval)))

    # flag entire flare region
    for ie in range(len(rss1[0])):
        flarefound[(xval >= (rss1[0])[ie]) & (xval <= ((rss1[1])[ie] + (rss1[0])[ie]))] = 1

    # if debug:
    #     di = [t - s for s, t in zip(xval, xval[1:])]
    #     xnogaps = [x - d + 0.1 for x,d in zip(xval[1:],di) if d > 0.5]
    #     plt.plot(xnogaps,yval,'k.')
    #     plt.plot(xnogaps[flarefound==1],yval[flarefound==1],'g.')
    #     plt.show()
    #     plt.close()

    allflaresfound = 0
    if np.size(rss1[0]) == 0:
        allflaresfound = 1

    run = 1
    while allflaresfound == 0:
        print("run again")
        run += 1

        lccut = deepcopy(yval)
        lccut[flarefound == 1] = replaceval

        lccutalt = deepcopy(yval)
        lccutalt = lccutalt[flarefound == 0]

        fluxdiffalt2 = (lccutalt - np.roll(lccutalt, 2))
        fluxdiffalt = (lccutalt - np.roll(lccutalt, 1))

        obs2 = xval[flarefound == 0]
        nd2 = daynight(obs2)

        # find new local standard deviation of fluxdiffs
        # p = Pool(12)

        geterrf = partial(geterr20, fluxdiffalt, nd2)
        # rmshp2 = p.map(geterrf, range(len(nd2) - 1))
        rmshp2 = []
        for r in range(len(nd2) - 1):
            rmshp2.append(geterrf(r))

        rmsh2 = np.hstack(np.array(rmshp2))
        if np.size(np.where(rmsh2 == 0)[0]) > 0:
            print("zero vals. check.")
            break

        # new fluxdiffs. flare-cut values replaced by median of fluxdiffs
        newfluxdiff = np.ones((np.size(yval)))
        newfluxdiff[flarefound == 0] = fluxdiffalt
        newfluxdiff[flarefound == 1] = np.nan#np.median(fluxdiffalt)

        newfluxdiff2 = np.ones((np.size(yval)))
        newfluxdiff2[flarefound == 0] = fluxdiffalt2
        newfluxdiff2[flarefound == 1] = np.nan#np.median(fluxdiffalt2)

        # new local rms of fluxdiffs. flare-cut rms values replaced by mean of rms
        rmsnew = np.ones((np.size(yval)))
        rmsnew[flarefound == 0] = rmsh2
        rmsnew[flarefound == 1] = np.nan#np.ma.mean(rmsh2)

        # find flares with new fluxdiffs and rms values.
        # lccut same as lc. but the flare-zone flux values are replaced by replaceval (0.9001). replaceval<1. -> ignored in the decay fit.
        rmsh_test = running_box(xval, lccut, 80 * exp, operation="std")

        rss8 = findflares_sub(thresh, sigma, lccut, xval, nd, newfluxdiff, newfluxdiff2, rmsnew,rmsh_test, exp, run,120,320,debug)

        # print(rss1)
        # print(rss8)

        # combine previously found flares with new flares
        rss1 = (np.hstack((rss1[0], rss8[0])), np.hstack((rss1[1], rss8[1])),np.hstack((rss1[2], rss8[2])),np.hstack((rss1[3], rss8[3])))
        # print(rss1)

        if np.size(rss8[0]) == 0:
            allflaresfound = 1
        else:
            i_peak = rss1[2]
            i_thresh = rss1[3]
            # print(i_peak,i_thresh)

        # print "size rss1: ",np.size(rss1[0])
        flarefound = np.zeros((np.size(xval)))

        for ie in range(len(rss1[0])):
            flarefound[(xval >= (rss1[0])[ie]) & (xval <= ((rss1[1])[ie] + (rss1[0])[ie]))] = 1

    return flarefound, i_peak, i_thresh


def load_test_npy():
    nr = 95
    xval = np.load("/appcg/data1/fl386/flaretest/xval_" + str(nr) + "_200120.npy")
    yval = np.load("/appcg/data1/fl386/flaretest/yval_" + str(nr) + "_200120.npy")
    return xval,yval

def load_lcs(lc_file):
    import csv
    jd, lc = [], []
    with open(lc_file, 'r') as sfile:
        reader = csv.reader(sfile)
        for row in reader:
            try:
                lc.append(float(row[1]))
                jd.append(float(row[0]))
            except Exception as e:
                print(e)
    sfile.close()
    lc = normalise(lc)
    print("EXPOSURE: " + str(24*60*60*(jd[1]-jd[0])) + " seconds")
    return np.array(jd), np.array(lc)

def normalise(a):
    return np.ma.divide(a, np.nanmean(a))

def main(x,y,thresh,sigma,debug=False):
    # used for plotting
    zeronr = np.int(np.min(x) / 1000.) * 1000.
    x -= zeronr
    nd = daynight(x)
    # just for plotting
    (t_axis, ndd, xticks) = getpl(x, nd, zeronr)

    flarefound, i_peak, i_thresh = findflares(x, y,thresh,sigma,debug)

    return flarefound, t_axis, i_peak, i_thresh


# if __name__ == "__main__":

    # # fname = "/appct/data/SPECULOOSPipeline/Io/output/lightcurves/Sp0440-0530/I+z/Sp0440-0530_6_v1globalunclip.csv"
    # # tel = "Callisto"
    # tel = "Io"
    # # targ = "Sp0721-3105"
    # targ = "Sp0440-0530"
    # ap = 6
    # thresh = 6
    # sigma = 3
    # plot = True
    # # targ = "Sp1626-3812"
    # fname = "/appct/data/SPECULOOSPipeline/" + tel + "/output/lightcurves/" + targ + "/I+z/" + targ + "_" + str(ap) +"_v1globalunclip.csv"
    # # x,y = load_test_npy()
    # x,y = load_lcs(fname)
    # # print(str(2400/(24*60*60*(jd[1]-jd[0]))))
    # flarefound, t_axis = main(x,y,thresh,sigma)
    #
    # if plot:
    #     plt.figure(figsize=(200, 5))
    #     plt.plot(t_axis[flarefound == 1], y[flarefound == 1], ".", color="red")
    #     plt.plot(t_axis[flarefound == 0], y[flarefound == 0], ".", color="green")
    #     # for i in i_peak:
    #     #     plt.axvline(x=t_axis[i],linestyle='--',color="k",alpha=0.5)
    #     #     print(t_axis[i])
    #
    #     plt.ylim(0.95, 1.05)#np.max(y))
    #     plt.savefig("/appct/data/SPECULOOSPipeline/tests/flares/florian/TEST_"+ targ + "_lowerthresh"+str(thresh)+"_flarezone120_sigma"+str(int(10*sigma))+"_j2_2exp")


