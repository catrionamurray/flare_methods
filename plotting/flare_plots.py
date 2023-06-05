from imports import *

def plot_binned_hist(x, y, z, xlabel, ylabel, zlabel, ytype, fname):
    from scipy.stats import binned_statistic_2d

    n = 10

    if ytype == 'log':
        y_bins = np.logspace(np.log10(np.nanmin(y)), np.log10(np.nanmax(y)), n)
    else:
        y_bins = np.linspace(np.nanmin(y), np.nanmax(y), n)

    x_bins = np.linspace(np.nanmin(x), np.nanmax(x), n)

    ret = binned_statistic_2d(x, y, z, statistic=np.mean, bins=[x_bins, y_bins])

    fig, (ax2) = plt.subplots(1, 1, figsize=(10, 8))
    im = ax2.pcolormesh(ret.x_edge, ret.y_edge, ret.statistic.T)
    cbar = fig.colorbar(im, ax=ax2)
    cbar.set_label(zlabel, rotation=270, labelpad=20)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel)
    if ytype == 'log':
        ax2.set_yscale('log')
    # plt.show()
    plt.savefig(fname)
    plt.close()


def plot_flare_dist(all_flares, flarelc, pltname):
    all_flares['energies'] = np.log10(flarelc.bol_lum) + 7 + all_flares['int']

    plt.scatter(all_flares['amp'], 24 * 60 * all_flares['fwhm'], c=all_flares['energies'])
    plt.xscale('log')
    # bins = [28, 28.5, 29, 29.5, 30, 30.5, 31,31.5,32,32.5,33,33.5,34,34.5]
    # all_flares['binned_energies'] = pd.cut(all_flares['energies'], bins=bins, labels=bins[:-1])
    #
    # for b in range(len(bins)-1):
    #     x=all_flares[all_flares['binned_energies'] == bins[b]]['amp']
    #     y=all_flares[all_flares['binned_energies'] == bins[b]]['fwhm']
    #     print(bins[b],len(x))
    #     plt.semilogx(x,y*24*60,'.',label=str(bins[b])+"-"+str(bins[b+1]))
    plt.xlabel("Amplitude")
    plt.ylabel("FWHM (minutes)")
    plt.axvline(x=0.01, linestyle='--', color='k')
    plt.xlim(0.00005, 10)
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.colorbar()
    # plt.show()
    plt.savefig(pltname + "_energies")
    plt.close()

    # plt.scatter(all_flares['amp'], 24*60*all_flares['fwhm'], c=all_flares['int'])
    # plt.xscale('log')
    # # print(max(all_flares['int']),min(all_flares['int']))
    # # bins = [-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6]
    # # all_flares['binned_int'] = pd.cut(all_flares['int'], bins=bins, labels=bins[:-1])
    # #
    # # for b in range(len(bins)-1):
    # #     x=all_flares[all_flares['binned_int'] == bins[b]]['amp']
    # #     y=all_flares[all_flares['binned_int'] == bins[b]]['fwhm']
    # #     plt.semilogx(x,y*24*60,'.',label=str(bins[b])+"-"+str(bins[b+1]))
    # plt.xlabel("Amplitude")
    # plt.ylabel("FWHM (minutes)")
    # plt.axvline(x=0.001,linestyle='--',color='r')
    # plt.xlim(0.00005, 10)
    # # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    # plt.colorbar()
    # # plt.show()
    # plt.savefig(pltname + "_int")
    # plt.close()


def plot_hist(x, bins, xlabel, ylabel, colour="Blue", xticks=None, svname=None):
    # plt.hist(x, bins=bins, edgecolor="white",facecolor=colour)
    x.plot.hist(bins=bins, color=colour)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    if xticks is not None: plt.xticks(xticks)
    if svname is None:
        plt.show()
    else:
        plt.savefig(svname)
    plt.close()


def plot_hist_rot(x1, x2, bins, xlabel, ylabel, label, svname=None, figsize=None):
    if figsize is not None:
        plt.figure(figsize=figsize)
    plt.hist(x1, bins=bins, facecolor="orange", edgecolor="white", label='')
    plt.hist(x2, bins=bins, facecolor="blue", edgecolor="white", label=label)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.legend(loc='upper right')
    if svname is None:
        plt.show()
    else:
        plt.savefig(svname)
    plt.close()


def plot_bar(x, y, xlabel, ylabel, xticks=None, width=1, colour="orange", svname=None):
    plt.bar(x=x, height=y, edgecolor="white", facecolor=colour, align='center', width=width)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    if xticks is not None: plt.xticks(xticks)
    if svname is None:
        plt.show()
    else:
        plt.savefig(svname)

    plt.close()

def get_fraction(x1, x2, bins):
    n1, b1, p1 = plt.hist(x1, bins=bins)
    n2, b2, p2 = plt.hist(x2, bins=bins)
    newb = [(a + b) / 2 for a, b in zip(b1[::1], b1[1::])]
    plt.close()

    return n2 / n1, newb
