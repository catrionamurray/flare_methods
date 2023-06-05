from imports import *
from utils import *
from statistics.flare_stats import calculate_bol_luminosity, flare_energy
from statistics.flare_inject_recovery import generate_flare_distribution, mod_random
from detection.flare_detection import flare_start_stop, get_flare_peak, get_flare_region, fit_davenport_model, qcheck_flares
from plotting.flare_plots import plot_flare_dist, plot_bar, plot_hist, plot_binned_hist, plot_hist_rot
from flare_rw import write_to_file, write_to_file_recovered, import_all_flares, create_flaredict
import detection.florian_flares as florian_flares
from modeling.aflare import aflare1

class FlareLightCurve:
    def __init__(self, t, f, teff, rad, id):
        self.t = t
        self.f = f
        self.r = rad
        self.teff = teff
        self.id = id

    def found_flares(self, flares_dict):
        """
        :param flares_dict:
        :return:
        """
        self.flares = pd.DataFrame(flares_dict)

    def mask_flares(self, flaremask):
        """
        Store the flux and time with flares masked
        :param flaremask:
        :return:
        """
        self.clean_f = self.f[flaremask == 0]
        self.clean_t = self.t[flaremask == 0]

    def calculate_bol_lum(self, instrument_response):
        # calculate flare energy
        response_fname = instrument_response  # "/appct/data/SPECULOOSPipeline/andorSPC_I+z_instrumentSR.csv"
        lam, response = get_response(response_fname)
        bol_lum = calculate_bol_luminosity(lam, response, self.teff, self.r)
        self.bol_lum = bol_lum

    def get_flares_energy(self,
                          amp_range,
                          fwhm_range,
                          mode_amp,
                          mode_fwhm,
                          energy_range,
                          e_bin,
                          num_flares=10,
                          outname=""):

        """
        Get the flare energies
        :param amp_range:
        :param fwhm_range:
        :param mode_amp:
        :param mode_fwhm:
        :param energy_range:
        :param e_bin:
        :param num_flares:
        :param outname:
        :return:
        """

        mintime, maxtime = np.min(self.t), np.max(self.t)
        tmodel = np.arange(min(self.t), max(self.t), 0.00001)

        # maximum and minimum flare limits
        flaremodel = aflare1(tmodel, np.median(tmodel), fwhm_range[0], amp_range[0])
        minenergy = np.log10(flare_energy(tmodel, flaremodel, self.bol_lum))
        flaremodel = aflare1(tmodel, np.median(tmodel), fwhm_range[1], amp_range[1])
        maxenergy = np.log10(flare_energy(tmodel, flaremodel, self.bol_lum))
        print("Minimum Flare Energy: %0.3f, Maximum Flare Energy: %0.3f " % (minenergy, maxenergy))

        flares_dict = {'peak_t': [], 'amp': [], 'fwhm': [], 'energy': []}

        energybins = np.arange(energy_range[0], energy_range[1], e_bin)
        energybins_range = [(energybins[i], energybins[i + 1]) for i in range(len(energybins) - 1)]
        print(energybins_range)

        for e in energybins_range:
            n = 0
            total = 0
            f_with_flare = self.f.copy()

            while n < num_flares:
                total = total + 1
                fwhm_dist, amp_dist = generate_flare_distribution(1, amp_range, fwhm_range, mode_amp, mode_fwhm)
                flaremodel = aflare1(tmodel, np.median(tmodel), fwhm_dist[0], amp_dist[0])
                energy = np.log10(flare_energy(tmodel, flaremodel, self.bol_lum))
                if energy >= e[0] and energy <= e[1]:
                    flares_dict['peak_t'].append(np.median(tmodel))
                    flares_dict['amp'].append(amp_dist[0])
                    flares_dict['fwhm'].append(fwhm_dist[0])
                    flares_dict['energy'].append(energy)
                    n = n + 1

            flaremodel = np.zeros(len(tmodel))

            for k in range(0, num_flares):
                ingap = True
                while ingap:
                    t_dist = (mod_random(1)[0] * (maxtime - mintime)) + mintime
                    nearest_t = find_nearest(self.t, t_dist)
                    # print(np.absolute(t_dist - self.t[nearest_t]),t_dist)
                    # ensure the surrounding +/- 2 points are within 0.01d (14.4 mins) of the flare to prevent flares in gaps/SON/EON
                    if (np.absolute(t_dist - self.t[nearest_t]) < 0.01) & (
                            np.absolute(t_dist - self.t[nearest_t - 2]) < 0.01) \
                            & (np.absolute(t_dist - self.t[nearest_t + 2]) < 0.01):
                        ingap = False

                tpeak = t_dist
                flares_dict['peak_t'][k] = tpeak
                flare = aflare1(self.t, tpeak, flares_dict['fwhm'][k], flares_dict['amp'][k])
                fmodel = aflare1(tmodel, tpeak, flares_dict['fwhm'][k], flares_dict['amp'][k])
                f_with_flare = f_with_flare + flare
                flaremodel = flaremodel + fmodel

            plt.figure(figsize=(24, 8))
            plt.title(e)
            plt.plot(tmodel, flaremodel, 'g.', alpha=0.01)
            plt.plot(self.t, f_with_flare - self.f, 'k.')
            plt.show()
            plt.close()

        self.f_inject = f_with_flare
        self.flares = pd.DataFrame(flares_dict)

    def inject_flare(self,
                     nfake,
                     repeat,
                     allflares,
                     energybins,
                     thresh,
                     sigma,
                     binsize,
                     outname,
                     pltname,
                     plot=True):
        """
        Inject synthetic flares into LCs
        :param nfake: number of synthetic flares to inject
        :param repeat: number of times to repeat the injection tests
        :param allflares:
        :param energybins:
        :param thresh:
        :param sigma:
        :param binsize: Binsize in minutes
        :param outname: File location for saving files
        :param pltname: File location for saving plots
        :return:
        """

        mintime, maxtime = np.min(self.t), np.max(self.t)
        tmodel = np.arange(min(self.t), max(self.t), 0.00001)
        targ_amps, targ_fwhms, targ_tpeaks, targ_energies, targ_recovered = [], [], [], [], []
        svname = pltname + "_injectedflares.csv"
        svname_all = outname + "allinjectedflares.csv"

        if os.path.exists(svname):
            os.remove(svname)
        if os.path.exists(svname.replace("injectedflares.csv", "recoveredflares.csv")):
            warnings.warn(f"Rerunning the injection-recovery tests will overwrite the {svname.replace('injectedflares.csv', 'recoveredflares.csv')} file!")
            os.remove(svname.replace("injectedflares.csv", "recoveredflares.csv"))

        # get RMS of each night, so that we can correlate with recovery rate
        stddev, bin_rms, error_within_bin = [], [], []
        ts, split = split_multilc(self.t) # split the global LC into separate nights
        fs = np.split(self.f, split)
        for t, f in zip(ts, fs): # loop over each night
            bin_t, bin_f, bin_e = bin_data(t, f, binsize)
            stddev.append(np.nanstd(f))
            bin_rms.append(np.nanstd(bin_f)) # get the standard deviation of binned flux
            error_within_bin.append(np.nanmean(bin_e)) # get the average error within each bin

        for en in range(len(energybins) - 1): # loop through each energy bin
            eflares = allflares[allflares['energy'] > energybins[en]]
            eflares = eflares[eflares['energy'] < energybins[en + 1]]
            eflares.reset_index(inplace=True, drop=True)
            print(str(len(eflares)), " flares in energy bin between ", str(energybins[en]), "and",
                  str(energybins[en + 1]))

            for r in range(repeat): # repeat the injection recovery tests for a user-defined number ("repeat") of times
                f_with_flares = self.f.copy()
                flare_fmodel = np.zeros(len(tmodel))

                if len(eflares['energy'].values) > nfake:
                    amps, fwhms, tpeaks, energies, rednoise, whitenoise = [], [], [], [], [], []

                    # loop over the number of fake flares you want to generate
                    for k in range(0, nfake):
                        ingap = True

                        x = mod_random(1)[0] * (len(eflares) - 1)
                        amp = eflares['amp'].values[int(round(x))]
                        fwhm = eflares['fwhm'].values[int(round(x))]

                        # if the injected flare occurs in a gap in data then generate a new flare.
                        while ingap:
                            t_dist = (mod_random(1)[0] * (maxtime - mintime)) + mintime
                            nearest_t = find_nearest(self.t, t_dist)
                            npoints = 5
                            # ensure the surrounding +/- 5 points are within 0.01d (14.4 mins) of the flare to
                            # prevent flares in gaps/SON/EON
                            if (nearest_t + npoints < len(self.t)) and (nearest_t - npoints >= 0):
                                if (np.absolute(t_dist - self.t[nearest_t]) < 0.01) & (
                                        np.absolute(t_dist - self.t[nearest_t - npoints]) < 0.01) \
                                        & (np.absolute(t_dist - self.t[nearest_t + npoints]) < 0.01):
                                    ingap = False

                        # generate synthetic flare from flare model
                        tpeak = t_dist
                        flare = aflare1(self.t, tpeak, fwhm, amp) # flare in LC
                        flaremodel = aflare1(tmodel, tpeak, fwhm, amp) # flare model for plotting/energy purposes

                        tpeaks.append(tpeak) # time of flare peak
                        amps.append(amp) # amplitude of flare
                        fwhms.append(fwhm) # fwhm of flare
                        energies.append(eflares['energy'].values[int(round(x))]) # flare energy

                        # loop over each day in the split LC
                        for ti in range(len(ts)):
                            # if the time of the flare peak is
                            if ts[ti][0] <= tpeak < ts[ti][-1]:
                                rednoise.append(bin_rms[ti])
                                whitenoise.append(error_within_bin[ti])

                        f_with_flares = f_with_flares + flare # add flare to LC
                        flare_fmodel = flare_fmodel + flaremodel # add flare model to model LC

                        eflares = eflares.drop(index=int(round(x)))
                        eflares.reset_index(inplace=True, drop=True)

                    targ_tpeaks.extend(tpeaks)
                    targ_amps.extend(amps)
                    targ_fwhms.extend(fwhms)
                    targ_energies.extend(energies)

                    if plot:
                        plt.figure(figsize=(24, 8))
                        plt.title(str(energybins[en]) + " - " + str(energybins[en + 1]))
                        for tp in tpeaks:
                            plt.axvline(tp, color='orange')
                        plt.plot(tmodel, flare_fmodel, 'g.', alpha=0.01)
                        plt.plot(self.t, f_with_flares - self.f, 'k.')
                        plt.savefig(pltname + "_flaremodel_energy" + str(energybins[en]) + "_n" + str(r))
                        plt.close()

                    self.f_inject = f_with_flares
                    # self.plot_split(self.f, self.f_inject,tpeaks,pltname + "_split_energy" + str(energybins[e])+ "_n" + str(r))

                    # FLARE DETECTION CODE ********************
                    flarefound, i_peak, i_thresh = self.flare_detection(thresh=thresh, sigma=sigma)
                    # flarefound, i_peak, i_thresh = self.flare_detection(thresh=5, sigma=2.5)
                    # flarefound2, i_peak2, i_thresh2 = self.flare_detection(thresh=6, sigma=3.0)
                    # flarefound3, i_peak3, i_thresh3 = self.flare_detection(thresh=6, sigma=2.5)
                    # flarefound4, i_peak4, i_thresh4 = self.flare_detection(thresh=5, sigma=3.0)
                    # *****************************************

                    i_peak = sorted(i_peak)
                    # get the start and stop times of flares
                    flares_start, flares_stop, flares_peak = flare_start_stop(flarefound,
                                                                              i_peak)  # get_flare_region(flarefound,i_peak)
                    flares_peak = get_flare_peak(flares_start, flares_stop, flares_peak, self.f_inject, self.t)

                    if len(i_peak) == 0:
                        print("Found no flares")
                        erec, trec, arec, frec, tpeaks_recover, amps_recover, fwhms_recover, energies_recover = [], [], [], [], [], [], [], []
                        recovered = np.zeros(nfake)
                        # all_recovered.extend(recovered)
                        targ_recovered.extend(recovered)
                        erec.extend(recovered)
                        trec.extend(recovered)
                        arec.extend(recovered)
                        frec.extend(recovered)

                    else:
                        # print("Recovered flares!")
                        # print("pause")

                        # for i in i_peak2:
                        #     plt.axvline(self.t[i]+0.001, color='k')
                        # for i in i_peak3:
                        #     plt.axvline(self.t[i] + 0.002, color='b')
                        # for i in i_peak4:
                        #     plt.axvline(self.t[i] + 0.003, color='orange')
                        # plt.show()
                        # plt.savefig(pltname + "_recovered_energy" + str(energybins[e]) + "_n" + str(r))
                        # plt.close()

                        recovered = []
                        # ts, split = split_multilc(self.t)
                        # fs = np.split(self.f_inject, split)

                        # x = find_nearest(ts, tp_inject)
                        fmask_temp = np.zeros(len(self.t))
                        energies_recover, tpeaks_recover, amps_recover, fwhms_recover, rms_day = [], [], [], [], []
                        # fmask = flarefound
                        for istart, istop in zip(flares_start, flares_start):
                            fmask_temp[istart:istart + 10] = 1

                        exp = self.t[1] - self.t[0]

                        for istart, istop, ipeak in zip(flares_start, flares_stop, flares_peak):
                            try:
                                # ok, tmodel2, fmodel2_median, fmodel2, popt,median_model,tflare,fflare = qcheck_flares(self.t, self.f_inject, istart,
                                #                                                         istop, ipeak, exp,
                                #                                                         fmask_temp,
                                #                                                         localsize=100,
                                #                                                         globalsize=320, npoints=2,
                                #                                                         nsigma=2.5, smooth=20,
                                #                                                         fittype="interp")

                                ok = True
                                popt = [np.nan, np.nan, np.nan]

                                if ok == True:
                                    if np.isnan(popt[0]):
                                        tpeaks_recover.append(self.t[ipeak])
                                    else:
                                        tpeaks_recover.append(popt[0])
                                    amps_recover.append(popt[2])
                                    fwhms_recover.append(popt[1])

                            except Exception as e:
                                print(e)
                                # tpeaks_recover.append(np.nan)
                                # amps_recover.append(np.nan)
                                # fwhms_recover.append(np.nan)
                                ok = False

                            if ok:
                                try:
                                    energy_rec = np.log10(
                                        flare_energy(tflare, (fflare / median_model) - 1,
                                                     self.bol_lum))

                                    # energy_rec = np.log10(
                                    #     flare_energy(self.t[istart:istop],self.f_inject[istart:istop]/median_model[istart:istop], self.bol_lum))
                                    energies_recover.append(energy_rec)
                                except Exception as e:
                                    print(e)
                                    energies_recover.append(np.nan)

                        # recovered_tpeak = []
                        if plot:
                            plt.figure(figsize=(24, 8))
                            plt.title(str(energybins[en]) + " - " + str(energybins[en + 1]))
                            plt.plot(tmodel, flare_fmodel, 'g.', alpha=0.01)
                            plt.plot(self.t, flare_f - self.f, 'k.')

                            for i in range(len(tpeaks)):
                                # j = i_peak[i]
                                # recovered_tpeak.append(self.t[j])
                                plt.axvline(tpeaks[i], color='r')
                                # plt.text(x=tpeaks[i] + 0.001, y=max(flare_fmodel), s=str(i_thresh[i]))

                        erec, trec, arec, frec = [], [], [], []

                        # for injected flare
                        for a_inject, f_inject, tp_inject in zip(amps, fwhms, tpeaks):
                            rec = False
                            # test if each recovered flare is within 30 minutes (to allow for multiple flares)
                            for tp_recover in range(len(tpeaks_recover)):
                                # print(np.absolute(t - tp))
                                # if np.absolute(tp_inject - tp_recover) < (10. / (24 * 60)):
                                # rec = True
                                if not np.isnan(tpeaks_recover[tp_recover]):
                                    if np.absolute(tp_inject - tpeaks_recover[tp_recover]) < f_inject:
                                        # print("RECOVERED WITHIN 1 FWHM")
                                        rec = True
                                        ind = tp_recover
                                    elif np.absolute(tp_inject - tpeaks_recover[tp_recover]) < (10. / (24 * 60)):
                                        print("more lenient condition")
                            if rec:
                                if plot:
                                    plt.axvline(tp_inject, color='g')
                                recovered.append(1.0)
                                targ_recovered.append(1.0)
                                erec.append(energies_recover[ind])
                                trec.append(tpeaks_recover[ind])
                                arec.append(amps_recover[ind])
                                frec.append(fwhms_recover[ind])
                                # else:
                                # rec= False
                                # recovered.append(0.0)
                                # targ_recovered.append(0.0)
                                # en = np.log10(integrate_flare(ts[x],fs[x]))
                            else:
                                recovered.append(0.0)
                                # all_recovered.append(0.0)
                                targ_recovered.append(0.0)
                                erec.append(0.0)
                                trec.append(0.0)
                                arec.append(0.0)
                                frec.append(0.0)

                        # plt.savefig(pltname + "_recovered_energy" + str(energybins[e]) + "_n" + str(r))
                        # plt.show()
                        if plot:
                            plt.savefig(pltname + "_recovered_energy" + str(energybins[en]) + "_n" + str(r))
                            plt.close()

                            self.plot_split(self.f, self.f_inject, tpeaks, recovered,
                                            pltname + "_split_energy" + str(energybins[en]) + "_n" + str(r))

                    if plot:
                        plot_binned_hist([24 * 60 * f for f in targ_fwhms], targ_amps, targ_recovered, "FWHM (mins)",
                                         "Amplitude", "Fraction Recovered", 'log', pltname + "_recovered_hist")
                        # plot_binned_hist([24 * 60 * f for f in all_fwhms], all_amps, all_recovered, "FWHM (mins)",
                        #                      "Amplitude", "Fraction Recovered",outname + "_recovered_hist")

                    write_to_file(tpeaks, amps, fwhms, energies, whitenoise, rednoise, trec, arec, frec, erec,
                                  recovered, svname)
                    write_to_file_recovered(tpeaks_recover, amps_recover, fwhms_recover, energies_recover,
                                            svname.replace("injectedflares.csv", "recoveredflares.csv"))
                    # write_to_file(tpeaks, amps, fwhms, energies, recovered, svname_all)

                else:
                    print("Less than 10 flares between " + str(energybins[en]) + " and " + str(energybins[en + 1]))

        return

    def inject_fake_flares(self, nfake, amp_range, fwhm_range, mode_amp, mode_fwhm):
        fwhm_dist, amp_dist = generate_flare_distribution(nfake, amp_range, fwhm_range, mode_amp, mode_fwhm)
        mintime, maxtime = np.min(self.t), np.max(self.t)
        flares_dict = {'peak_t': [], 'amp': [], 'fwhm': [], 'energy': []}
        flare_f = self.f.copy()
        tmodel = np.arange(min(self.t), max(self.t), 0.00001)
        flare_fmodel = np.zeros(len(tmodel))
        flare_fs, flare_ts, flare_fsbase = [], [], []

        # loop over the numer of fake flares you want to generate
        for k in range(0, nfake):
            ingap = True

            while ingap:
                t_dist = (mod_random(1)[0] * (maxtime - mintime)) + mintime
                nearest_t = find_nearest(self.t, t_dist)
                # ensure the surrounding +/- 2 points are within 0.01d (14.4 mins) of the flare to prevent flares in gaps/SON/EON
                if (np.absolute(t_dist - self.t[nearest_t]) < 0.01) & (
                        np.absolute(t_dist - self.t[nearest_t - 2]) < 0.01) \
                        & (np.absolute(t_dist - self.t[nearest_t + 2]) < 0.01):
                    ingap = False

            tpeak = t_dist
            flare = aflare1(self.t, tpeak, fwhm_dist[k], amp_dist[k])
            flaremodel = aflare1(tmodel, tpeak, fwhm_dist[k], amp_dist[k])

            energy = flare_energy(self.t, flare, self.bol_lum)
            energymodel = flare_energy(tmodel, flaremodel, self.bol_lum)
            print(np.log10(energy), np.log10(energymodel))
            flare_f = flare_f + flare
            flare_fmodel = flare_fmodel + flaremodel

            flares_dict['peak_t'].append(t_dist)
            flares_dict['amp'].append(amp_dist[k])
            flares_dict['fwhm'].append(fwhm_dist[k])
            flares_dict['energy'].append(np.log10(energymodel))

            # flare_fs.extend(flare_f[flare_fistart:istop])
            # flare_fsbase.extend(flare_f[istart:istop]-self.f[istart:istop])
            # flare_ts.extend(self.t[istart:istop])

        plt.figure(figsize=(24, 8))
        plt.plot(self.t, self.f, 'k.')
        plt.plot(self.t, flare_f, 'r.', alpha=0.35)
        plt.ylim(0.98, 1.02)
        plt.show()
        plt.close()

        plt.figure(figsize=(24, 8))
        plt.plot(tmodel, flare_fmodel, 'g.', alpha=0.01)
        plt.plot(self.t, flare_f - self.f, 'k.')
        plt.show()
        plt.close()

        self.f_inject = flare_f
        self.flares = pd.DataFrame(flares_dict)

    def flare_detection(self, thresh, sigma):
        flux = np.ma.array(self.f_inject.copy())
        time = self.t.copy()
        # with warnings.catch_warnings():
        #     warnings.simplefilter("ignore")
        flarefound, t_axis, i_peak, i_thresh = florian_flares.main(time, flux, thresh, sigma, debug=False)

        print("flare recovery done")
        return flarefound, i_peak, i_thresh

    def plot_split(self, y1, y2, tpeaks, rec, pltname):
        ts, split = split_multilc(self.t)
        y2s = np.split(y2, split)
        # ts, split = split_multilc(self.t)
        y1s = np.split(y1, split)
        # *********************
        fig, ax = plt.subplots(ncols=len(split) + 1, nrows=2, figsize=(42, 8), sharey=True, sharex=True)
        plt.subplots_adjust(wspace=0.0)
        if len(split) > 1:
            for i in range(len(split) + 1):
                for tp in range(len(tpeaks)):
                    if tpeaks[tp] >= ts[i][0] and tpeaks[tp] < ts[i][-1]:
                        if rec[tp] == 1.0:
                            ax[1][i].axvline(math.modf(tpeaks[tp])[0], color='g')
                        else:
                            ax[1][i].axvline(math.modf(tpeaks[tp])[0], color='r')
                tminus = [math.modf(j)[0] for j in ts[i]]
                ax[0][i].plot(tminus, y1s[i], 'k.')
                ax[1][i].plot(tminus, y2s[i], 'k.')
        else:
            for tp in tpeaks:
                ax[1][0].axvline(math.modf(tp)[0], color='r')
            tminus = [math.modf(j)[0] for j in ts[0]]
            ax[0][0].plot(tminus, y1s[0], 'k.')
            ax[1][0].plot(tminus, y2s[0], 'k.')
        # *********************
        plt.tight_layout()
        # plt.show()
        plt.savefig(pltname)
        plt.close()

