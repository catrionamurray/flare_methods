from imports import *
from flare_class import FlareLightCurve
from flare_rw import import_all_flares
from plotting.flare_plots import plot_flare_dist


def mod_random(x, d=False, seed=667):
    """
    Helper function that generates deterministic
    random numbers if needed for testing.
    Parameters
    -----------
    d : False or bool
        Flag to set if random numbers shall be deterministic.
    seed : 5 or int
        Sets the seed value for random number generator.
    """
    if d == True:
        np.random.seed(seed)
        return np.random.rand(x)
    else:
        np.random.seed()  # do not remove: seed is fixed otherwise!
        return np.random.rand(x)

def generate_all_flares(t_range, amp_range, fwhm_range, numflares, outname, plot=True):
    # tmodel = np.arange(t_range[0], t_range[1], 0.00001)
    all_flares = {'amp': [], 'fwhm': [], 'int': []}
    fwhm_dist, amp_dist = generate_flare_distribution(numflares, amp_range, fwhm_range, "lognormal", 'uniform')

    bins = np.arange(amp_range[0], 0.1, 0.001)
    if plot:
        plt.hist(amp_dist, bins=bins, histtype='stepfilled', alpha=0.4, label="lognormal")
        plt.legend(loc='best', frameon=False)
        plt.savefig(outname.replace(".csv", "_amp.png"))
        plt.close()
    # count = 0

    # en = []
    # for f, a in zip(fwhm_dist, amp_dist):
    #     if count%10 ==0:
    #         print(count)
    #     count = count + 1
    #     flaremodel = aflare1(tmodel, np.median(tmodel), f, a)
    #     # energy = np.log10(flare_energy(tmodel, flaremodel, self.bol_lum))
    #     en = np.log10(integrate_flare(tmodel, flaremodel))  # + np.log10(self.bol_lum) + 7
    #     all_flares['amp'].append(a)
    #     all_flares['fwhm'].append(f)
    #     all_flares['int'].append(en)
    #
    # plt.hist(all_flares['int'], histtype='stepfilled', alpha=0.4, label="flare energies")
    # plt.legend(loc='best', frameon=False)
    # plt.savefig(outname.replace(".csv","_int.png"))
    # plt.close()

    if outname != "":
        all_flares = pd.DataFrame(all_flares)
        all_flares.to_csv(outname, index=False)
        print(outname)

    return all_flares


def generate_flare_distribution(nfake, amp=[1e-4, 5], fwhm=[0.0025, 0.012], mode_a="uniform", mode_f="uniform",
                                **kwargs):
    """
    ADAPTED FROM ALTAIPONY
    Creates different distributions of fake flares to be injected into light curves.
    "uniform": Flares are distibuted evenly in duration and amplitude space.

    Parameters
    -----------
    nfake: int
        Number of fake flares to be created.
    ampl: [1e-4, 1e2] or list of floats
        Amplitude range in relative flux units.
    fwhm: [10, 2e4] or list of floats
        Duration range in days.
    mode: 'uniform'
        Distribution of fake flares in (duration, amplitude) space.
    kwargs : dict
        Keyword arguments to pass to mod_random
    Return
    -------
    dur_fake: durations of generated fake flares in days
    ampl_fake: amplitudes of generated fake flares in relative flux units
    """

    def generate_range(n, tup, **kwargs):
        x = (mod_random(n, **kwargs) * (tup[1] - tup[0])) + tup[0]
        return x

    def generate_normal_range(n, tup, **kwargs):
        x = np.random.normal(tup[0], (tup[1] - tup[0]) / 3., 1000)

        for i in range(n):
            j = x[i]
            while j < tup[0] or j > tup[1]:
                j = np.random.normal(tup[0], (tup[1] - tup[0]) / 2., size=1)
            x[i] = j

        return x

    def generate_lognormal_range(n, tup, **kwargs):
        sigma = 15
        mean = 0
        x = (np.random.lognormal(mean=mean, sigma=sigma, size=n))

        for i in range(n):
            j = x[i]
            while j < tup[0] or j > tup[1]:
                j = np.random.lognormal(mean=mean, sigma=sigma, size=1)
            x[i] = j

        return x

    if mode_f == 'uniform':
        fwhm_fake = generate_range(nfake, fwhm, **kwargs)

    if mode_a == 'uniform':
        ampl_fake = generate_range(nfake, amp, **kwargs)

    if mode_a == 'lognormal':
        ampl_fake = generate_lognormal_range(nfake, amp, **kwargs)

    if mode_a == 'normal':
        ampl_fake = generate_normal_range(nfake, amp, **kwargs)

    return fwhm_fake, ampl_fake

def inject_recovery_target(t, f, teff, rad, target, C, max_num=5, repeat=5, thresh=6, sigma=3.0, outname="",
                           pltname=""):
    # numflares =mod_random(max_num)[0]

    flarelc = FlareLightCurve(t, f, teff, rad, target)

    # generate a pool of random flare amps and fwhms
    # all_flares = generate_all_flares([min(flarelc.t),max(flarelc.t)],[0.0001,5.0],[0.5/(24*60),1./24],100000,os.path.dirname(outname) + "/all_flares_100000.csv")
    # quit()
    # import random flares
    all_flares = import_all_flares(os.path.dirname(outname) + "/all_flares_100000.csv")

    # calculate the bolometric luminosities of each star, dependant on radius and teff
    # flarelc.calculate_bol_lum()
    flarelc.bol_lum = 10 ** C  # 5.679776688240452e+22
    print(np.log10(flarelc.bol_lum) + 7)
    print("Thresh ", thresh, "Sigma", sigma, "Repeat ", repeat, " Max Num ", max_num)

    # all_flares = all_flares[24*all_flares['fwhm']<=1]

    plot_flare_dist(all_flares, flarelc, pltname)

    # inject fake flares into LC
    all_flares['energy'] = all_flares['int'] + C + 7
    all_flares = all_flares.sort_values(by='energy')
    # flarelc.inject_fake_flares(max_num,amp_range=[0.001,1.0],fwhm_range=[0.00139,0.012],mode_amp='lognormal',mode_fwhm='uniform')
    # flarelc.get_flares_energy(amp_range=[0.0001,1.0], fwhm_range=[1./(24*60),2./24], mode_amp='lognormal', mode_fwhm='uniform', energy_range=[29.0,34.0], e_bin=0.5, num_flares=max_num,outname=outname)
    # energybins = np.arange(28,35,1)
    energybins = np.arange(28, 35, 1)
    flarelc.inject_flare(max_num, repeat, all_flares, energybins, thresh, sigma, outname, pltname)

    return flarelc