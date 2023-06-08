from imports import *
import celerite
from celerite import terms
from scipy.optimize import minimize
from utils import *


def clip_data_for_gp(t,f,e,binsize,running_mean_box):
    bin_t, bin_f, bin_e = bin_data(t, f, binsize, e)
    nsigma_upper = 3
    nsigma_lower = 3
    bin_f = sigma_clip(bin_f, sigma_upper=nsigma_upper, sigma_lower=nsigma_lower).filled(np.nan)

    # print("Nsigma = ", nsigma_lower * np.nanstd(bin_f))
    # bin_f = sigma_clip(bin_f, sigma_upper=nsigma_upper, sigma_lower=nsigma_lower).filled(np.nan)

    cm = len(bin_f)  # 0
    it = True
    bin_f_clip = bin_f.copy()
    bin_t_clip = bin_t.copy()
    bin_e_clip = bin_e.copy()
    while it == True:
        run_med = running_box(bin_t_clip, bin_f_clip, 2*running_mean_box, 'median')
        run_std = running_box(bin_t_clip, bin_f_clip, 2*running_mean_box, 'std')
        avg_std = np.nanmedian(run_std)
        cond = np.logical_and(bin_f_clip < run_med + (nsigma_upper * avg_std),
                              bin_f_clip > run_med - (nsigma_lower * avg_std))
        bin_e_clip = bin_e_clip[cond]
        bin_t_clip = bin_t_clip[cond]
        bin_f_clip = bin_f_clip[cond]
        cmasked = len(bin_f_clip)
        if cmasked - cm == 0:
            it = False
        else:
            cm = cmasked

    bin_f = bin_f_clip
    bin_t = bin_t_clip
    bin_e = bin_e_clip

    # norm_t, norm_f = format_xy_plot(globallc.clean_t.copy(), globallc.clean_f.copy())
    #
    # yerr = bin_e[np.logical_and(~np.isnan(bin_f), ~np.isnan(bin_e))]
    # x = bin_t[np.logical_and(~np.isnan(bin_f), ~np.isnan(bin_e))]
    # y = bin_f[np.logical_and(~np.isnan(bin_f), ~np.isnan(bin_e))]

#     gp_f, gp_pred, wn = gp(x, y - 1, sq_err(yerr), x_pred, globallc, plot=False, figsize=(16, 4))
#                         #  gp_f2,gp_f3,gp_f4
#                         # gp_f2, gp_pred2 = gp(x, y - 1, mult_err(yerr,1.5), x_pred, globallc, plot=False, figsize=(16, 4))
#                         # gp_f3, gp_pred3 = gp(x, y - 1, add_err(yerr,0.02), x_pred, globallc, plot=False, figsize=(16, 4))
#                         # gp_f4, gp_pred4 = gp(x, y - 1, sq_err(yerr), x_pred, globallc, plot=False, figsize=(16, 4))
#
#                         # norm_bin_t, norm_bin_f = format_xy_plot(bin_t.copy(), bin_f.copy())
#                         norm_bin_t, norm_bin_f, norm_bin_e = bin_data(norm_t, globallc.clean_f, binsize,
#                                                                       globallc.clean_e)
#
#                         plt.figure(figsize=(16, 4))
#                         plt.plot(norm_t, globallc.clean_f, marker='.', color='c', linestyle="None")
#                         plt.plot(norm_bin_t, norm_bin_f, marker='.', color='k', linestyle="None")
#                         plt.title("Entire LC for " + t)
#                         plt.ylim(0.97, 1.03)
#                         if show_plots:
#                             plt.show()
#                         else:
#                             plt.savefig(pltname + '_global')
#                         plt.close()
#
#                         plt.figure(figsize=(16, 4))
#                         plt.plot(norm_t, globallc.clean_f, marker='.', color='c', linestyle="None")
#                         plt.plot(norm_bin_t, norm_bin_f, marker='.', color='k', linestyle="None")
#                         plt.plot(norm_t, gp_f, color='orange', lw=2)
#                         # plt.plot(norm_t, gp_f2, color='green',lw=2)
#                         # plt.plot(norm_t, gp_f3, color='red',lw=2)
#                         # plt.plot(norm_t, gp_f4, color='blue', lw=2)
#                         plt.title("Entire (with GP) LC for " + t)
#                         plt.ylim(0.97, 1.03)
#                         if show_plots:
#                             plt.show()
#                         else:
#                             plt.savefig(pltname + '_gp')
#                         plt.close()
#
#                         plt.figure(figsize=(24, 4))
#                         plt.plot(globallc.clean_t, globallc.clean_f, marker='.', color='c', linestyle="None")
#                         plt.plot(bin_t, bin_f, marker='.', color='k', linestyle="None")
#                         plt.plot(globallc.clean_t, gp_f, color='orange')
#                         plt.plot(x_pred, gp_pred, color='orange')
#                         plt.title("Entire (with GP) LC for " + t)
#                         plt.ylim(0.97, 1.03)
#                         if show_plots:
#                             plt.show()
#                         else:
#                             plt.savefig(pltname + '_gp_t')
#                         plt.close()
#
#                         if show_plots:
#                             norm_t = lc_plot(x=globallc.clean_t, y=globallc.clean_f / gp_f, dy=globallc.clean_e,
#                                              binsize=binsize, title="Entire (GP-detrended) Global LC for " + t)
#                         else:
#                             norm_t = lc_plot(x=globallc.clean_t, y=globallc.clean_f / gp_f, dy=globallc.clean_e,
#                                              binsize=binsize, title="Entire (GP-detrended) Global LC for " + t,
#                                              svname=pltname + '_global_gpdetrend')

def nllGP(p,y,t,e,gp):
    # if p[-1] < -1:
    #     return 1e25
    # print(p)
    gp.set_parameter_vector(p)
    try:
        gp.compute(t, yerr=e)
    except:
        return 1e25
    return -gp.log_likelihood(y)

def get_GP(gpk,kernel,t,y,e,am):
    # am_model = AirmassModel(a=-0.1,b=1)

    gp = celerite.GP(kernel, mean=np.nanmean(y), fit_mean=False)
    nmax=100000
    n = len(t)
    if n>nmax:
        l = np.sort(np.random.choice(n,nmax,replace=False))
    else:
        l = np.arange(n).astype(int)

    p0 = gp.get_parameter_vector()
    bounds = gp.get_parameter_bounds()
    # print('Initial GP HPs:', p0)
    gp.compute(t, yerr=e)
    # print('Initial NLL:', -gp.lnlikelihood(y))
    # soln = minimize(nllGP, p0)
    soln = minimize(nllGP, p0, method="L-BFGS-B", bounds=bounds, args=(y,t,e,gp))
    p1 = soln.x
    print('Initial GP HPs:', p0,'Fitted GP HPs:', p1)
    if gpk == "SHO":
        period = 2*math.pi/math.exp(p1[1])
        print('\nFitted Params: S0 = ',str(math.exp(p1[0])),'\nw0 = ',str(math.exp(p1[1])),"\nPeriod = ",period," days\n")
    elif gpk == "M32":
        print('\nFitted Params: Amp = ', str(np.sqrt(np.exp(p1[1]))), ", ls = " , str(np.exp(p1[0])))
    gp.set_parameter_vector(p1)
    gp.compute(t, yerr=e)
    mu = gp.predict(y, t, return_var=False, return_cov=False)
    r = y / mu

    return r, mu

def get_kernel(gpkernel,w0,sigma):

    if gpkernel == "SHO":
        Q = 1.0 / np.sqrt(2.0)
        # w0 = 3.0
        S0 = sigma ** 2 / (w0 * Q)
        print(np.log(w0), np.log(S0), np.log(Q))

        bounds = dict(log_S0=(-30, 30), log_Q=(-30, 30), log_omega0=(-30, np.log(2 * const.pi / 0.1)))

        kernel = terms.SHOTerm(log_S0=np.log(S0), log_Q=np.log(Q), log_omega0=np.log(w0),
                                bounds=bounds)
        kernel.freeze_parameter("log_Q")  # We don't want to fit for "Q" in this term

    elif gpkernel == "M32":
        # sigma = np.nanstd(y)#**2
        rho = 0.5/w0
        print(rho, np.log(rho))
        print(sigma, np.log(sigma))
        print(w0)
        print(0.5/w0)

        bounds = dict(log_rho=(-30,30),log_sigma=(-30,30))

        kernel = terms.Matern32Term(log_rho=np.log(rho),log_sigma=np.log(sigma),bounds=bounds) # * np.var(y)
        # kernel.freeze_parameter("log_sigma")
        # print('GP par: Amp = {0:.5f}, ls = {1:.5f}'.format(np.sqrt(np.exp(GP_par[0])), np.exp(GP_par[1])))
    return kernel


def gauss_proc(gpkernel,t,y,dy,airmass,added_sigma,w0,pltname):
    # added_sigma = 0.007
    added_sigma2 = 0.015
    added_sigma3 = 0
    # bounds = dict(log_S0=(-30, 30), log_Q=(-30, 30), log_omega0=(-30, np.log(2*const.pi/0.1)))

    try:
        print(np.var(y), np.nanstd(y)**2)
    except:
        pass

    kernel = get_kernel(gpkernel,w0,np.nanstd(y))

    sigma = [math.sqrt(i**2 + added_sigma**2) for i in dy]
    sigma2 = [math.sqrt(i**2 + added_sigma2**2) for i in dy]
    sigma3 = [math.sqrt(i ** 2 + added_sigma3 ** 2) for i in dy]

    # gp = celerite.GP(kernel, mean=np.nanmedian(y))
    # gp.compute(t, sigma)  # You always need to call compute once.

    r,mu = get_GP(gpkernel, kernel,t,y, sigma, airmass)
    r2, mu2 = get_GP(gpkernel, kernel, t, y, sigma2, airmass)
    r3, mu3 = get_GP(gpkernel, kernel, t, y, sigma3, airmass)

    # t,f = format_xy_plot(t.copy(),y)
    # x,pred_mean_f = format_xy_plot(x,pred_mean)

    color = "#ff7f0e"
    plt.figure(figsize=(200, 5))
    plt.plot(t, y, "k.", alpha=0.3,zorder=3)
    bin_t,bin_y,bin_e = bin_data(t,y,10,dy)
    plt.plot(bin_t,bin_y, "k.",zorder=7)
    # plt.errorbar(t, y, yerr=sigma2, fmt=".",color='g', capsize=0, alpha=0.25,zorder=0)
    # plt.errorbar(t, y, yerr=sigma, fmt=".", color=color,capsize=0, alpha=0.25,zorder=1)
    # plt.errorbar(t, y, yerr=sigma3, fmt=".",color='c', capsize=0, alpha=0.25,zorder=2)
    tdiff = np.positive(t-np.roll(t, 1))
    # print(tdiff)
    # print(max(tdiff),min(tdiff))
    plt.fill_between(t, y + sigma2, y - sigma2, where=tdiff<0.001, color='yellow', alpha=0.2,
                     edgecolor="none",zorder=0)
    plt.fill_between(t, y + sigma, y - sigma,where=tdiff<0.001, color=color,alpha=0.25,
                     edgecolor="none",zorder=1)
    plt.fill_between(t, y + sigma3, y - sigma3,where=tdiff<0.001, color='aqua',alpha=0.4,
                     edgecolor="none",zorder=2)

    plt.plot(t, mu3, color='aqua', label="Added sigma = " + str(added_sigma3),alpha=0.8,zorder=4)
    plt.plot(t, mu, color=color,label="Added sigma = " + str(added_sigma),zorder=6)
    plt.plot(t, mu2, color='yellow', label="Added sigma = " + str(added_sigma2),alpha=0.8,zorder=5)
    plt.xlabel("BJD")
    plt.ylabel("Relative Flux")
    plt.ylim(0.95,1.05)
    plt.legend()
    plt.savefig(pltname+"_gp")
    plt.close()

    # plot detrended version
    detrend_y = r#y/pred_mean
    plt.figure(figsize=(200, 5))
    plt.plot(t, detrend_y, "k.")
    # plt.errorbar(t, detrend_y, yerr=dy, fmt=".k", capsize=0, alpha=0.3)
    # plt.plot(t, pred_mean, color=color)
    plt.xlabel("BJD")
    plt.ylabel("Relative Flux")
    plt.ylim(0.95, 1.05)
    plt.savefig(pltname+"_gp_detrend")
    plt.close()

    return detrend_y