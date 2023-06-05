from imports import *
from utils import *

def flare_amp(f_array, t_array, t):
    ind = find_nearest(t_array, t)
    return f_array[ind]


def integrate_flare(t, f):
    return integrate.simps(f, t * 24 * 3600)


def flare_energy(t, f, bol_luminosity):
    energy = bol_luminosity * integrate_flare(t, f)  # W s = J
    energy_ergs = energy / 0.0000001
    return energy_ergs

def calculate_bol_luminosity(lam, response, teff, rad):
    teff_flare = 9000  # K
    # teff = 2700 # K

    r_star = rad * 695700 * 1000  # m

    planck_star = [pyasl.planck(teff, x) for x in lam]  # W/(m^2 m)
    planck_flare = [pyasl.planck(teff_flare, x) for x in lam]  # W/(m^2 m)
    int_star = integrate.simps(response * planck_star, lam)
    int_flare = integrate.simps(response * planck_flare, lam)

    A_flare = const.pi * (r_star ** 2) * int_star / int_flare  # m^2
    bol_luminosity = const.Stefan_Boltzmann * (teff_flare ** 4) * A_flare  # W

    return bol_luminosity


def abiogenesis_zone(E_u, R, T):
    flare_freq = [25.5 * ((10 ** 34) / (10 ** E)) * (R ** 2) * (T ** 4) for E in E_u]
    return flare_freq


def energy_uband(int_flare, teff_flare, bolE):
    # teff_flare = 9000  # K
    # planck_flare = [pyasl.planck(teff_flare, x) for x in lam]  # W/(m^2 m)
    # int_flare = integrate.simps(response * planck_flare, lam)
    A_flare = (10 ** bolE) / (const.Stefan_Boltzmann * (teff_flare ** 4))
    E_u = A_flare * int_flare
    print(int_flare)
    return np.log10(E_u)


def ffd(uniquetargs, targets, energy, rad, teff, total_t, int_flare, teff_flare, numflares=2, svname=None):
    thresh_e = []
    thresh_inds = []
    plt.figure()
    for t in range(len(uniquetargs)):
        inds = np.where(np.array(targets) == uniquetargs[t])[0]
        #     print(t, inds)
        energies = np.array(energy)[inds]
        #     print(energies)
        energies = [x for x in energies if ~np.isnan(x)]
        if len(energies) >= numflares:
            thresh_e.append(energies)
            thresh_inds.append(inds)

            logEs = np.arange(min(energies) - 0.1, max(energies) + 0.1, 0.1)

            logE_us = []
            for E in logEs:
                logE_u = 0.12666 * E  # energy_uband(int_flare,teff_flare, E)
                logE_us.append(logE_u)

            # logE_us = np.arange(energy_uband(response_u, lam_u, min(energies) - 0.1), energy_uband(response_u, lam_u, max(energies) + 0.1), 0.1)
            flarerate = cumulative_flare_rate(energies, logEs, total_t[inds[0]])
            ab_flarerate = abiogenesis_zone(logE_us, rad[inds[0]], teff[inds[0]] / 5778.)

            x = logEs[flarerate != 0.0]
            y = flarerate[flarerate != 0.0]
            points = plt.semilogy(x, y, marker=".", linestyle="None")
            abzone = plt.semilogy(logEs, ab_flarerate, color="green")
            col = points[0].get_color()
            m, b = np.polyfit(x, np.log10(y), 1)
            x_line = np.linspace(28, 36, 1000)
            y_line = m * x_line + b
            plt.semilogy(x_line, [np.power(10, yi) for yi in y_line], color=col, alpha=0.3)
            print("alpha: " + str(m), "beta: " + str(b))

    plt.xlim(28, 36)
    plt.ylim(0.001, 20)
    plt.ylabel("Flare Rate > logE (d^-1)")
    plt.xlabel("logE")
    if svname is None:
        plt.show()
    else:
        plt.savefig(svname)

    plt.close()

    return thresh_e, thresh_inds

def flare_rate(t, num_flares):
    total_time = get_total_observed_time(t)
    per_day = num_flares / total_time
    return per_day

def num_flare_type(flaretype, ftype="F"):
    return np.where(np.array(flaretype) == ftype)[0]

def get_flare_rate(t, num_flares):
    return num_flares / t

def get_num_flares(targets):
    targs = Counter(targets).keys()  # equals to list(set(words))
    numflares = Counter(targets).values()  # counts the elements' frequency

    return targs, numflares

def cumulative_flare_rate(flare_logEs, logEs, time):
    #     print(flare_logEs,time)
    # logEs = range(min(flare_logEs),max(flare_logEs),step=0.1)

    sorted_flare_logEs = sorted(flare_logEs)
    # get cumulative frequency of the flare energies
    res = stats.cumfreq(sorted_flare_logEs, numbins=len(logEs), defaultreallimits=(min(logEs), max(logEs)),
                        weights=None)
    #     print(res.cumcount)
    #     print(max(res.cumcount)-res.cumcount)
    # print(res.binsize)
    freq = (max(res.cumcount) - res.cumcount) / time

    return freq


