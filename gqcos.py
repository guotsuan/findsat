"""
This is cosmology module written by GQ

>>> abs_m(17.7, 0.1)
array([-20.61520457])
"""

import numpy as np
import scipy.stats as stat
import cosmolopy.distance as cd
import cosmolopy.constants as cc
from scipy.interpolate import interp1d
from scipy.stats import gaussian_kde
from sklearn.neighbors import KernelDensity
from statsmodels.nonparametric.kde import KDEUnivariate
from scipy.integrate import quad


# cosmo = {'omega_M_0': 0.3, 'omega_lambda_0': 0.7, 'omega_k_0': 0.0, 'h':
# 0.70}

def bootstrap(sample, ns=None, return_index=False):
    """
    Arguments:
       sample - input sample of values
       ns - number of samples to generate

    Performs resampling from sample with replacement, gathers
    statistic in a list computed by statfunc on the each generated sample.
    """
    if ns is None:
        ns = np.asarray(sample).shape[0]
    # print "input sample = ",  ns
    nsample = np.asarray(sample).shape[0]
    resample = [sample[j] for j in stat.randint.rvs(0, nsample - 1, size=ns)]

    if return_index:
        return [j for j in stat.randint.rvs(0, nsample - 1, size=ns)]
    else:
        return resample


def vr(halom):
    h0 = 0.72
    delta = 200
    factor = 4.0 * np.pi / 3.0
    rho_mean = 2.7752e11 * h0 * h0
    return 1000 * (halom / factor / delta / rho_mean) ** (1.0 / 3.0)


def sm_00(g, r):
    r0 = 4.64
    lr = 10 ** (-0.4 * (r - r0))
    x = g - r
    # a = 2.277
    # b = -0.398794
    a = 3.11971
    b = -0.546293
    result = (a * x + b) * lr
    return result


def abs_m(m, z, **cosmo):
    """ Return the absolute magnitude
    m: apparent magnitude
    z: redshift

    >>> abs_m(17.7, 0.1)
    array([-20.61520457])
    """

    m = np.atleast_1d(m)
    z = np.atleast_1d(z)
    if cosmo == {}:
        cosmo_in = {'omega_M_0': 0.3, 'omega_lambda_0': 0.7,
                    'omega_k_0': 0.0, 'h': 0.70}
    else:
        cosmo_in = cosmo

    lum_dist = cd.luminosity_distance(z, **cosmo_in)

    absm = - 25. + m - 5.0 * np.log10(lum_dist)
    return absm


def abs_mf(m, z, **cosmo):
    m = np.atleast_1d(m)
    z = np.atleast_1d(z)
    if cosmo == {}:
        cosmo_in = {'omega_M_0': 0.3, 'omega_lambda_0': 0.7,
                    'omega_k_0': 0.0, 'h': 0.70}
    else:
        cosmo_in = cosmo

    lum_dist = cd.quick_distance_function(cd.luminosity_distance, **cosmo_in)
    absm = - 25. + m - 5.0 * np.log10(lum_dist(z))
    return absm


def lum_distf(z, **cosmo):
    z = np.atleast_1d(z)
    if cosmo == {}:
        cosmo_in = {'omega_M_0': 0.3, 'omega_lambda_0': 0.7,
                    'omega_k_0': 0.0, 'h': 0.70}
    else:
        cosmo_in = cosmo

    print cosmo_in
    lum_dist = cd.quick_distance_function(cd.luminosity_distance, **cosmo_in)

    return lum_dist(z)


def lb_time(z, **cosmo):
    z = np.atleast_1d(z)
    if cosmo == {}:
        cosmo_in = {'omega_M_0': 0.3, 'omega_lambda_0': 0.7,
                    'omega_k_0': 0.0, 'h': 0.70}
    else:
        cosmo_in = cosmo

    results = cd.lookback_time(z, **cosmo_in) / 60./60./24./365./1.e9
    return results


def age_t(z, **cosmo):
    z = np.atleast_1d(z)
    if cosmo == {}:
        cosmo_in = {'omega_M_0': 0.3, 'omega_lambda_0': 0.7,
                    'omega_k_0': 0.0, 'h': 0.70}
    else:
        cosmo_in = cosmo

    results = cd.age(z, **cosmo_in) / cc.Gyr_s
    return results


def comoving_dist(z, **cosmo):
    z = np.atleast_1d(z)
    if cosmo == {}:
        cosmo_in = {'omega_M_0': 0.3, 'omega_lambda_0': 0.7,
                    'omega_k_0': 0.0, 'h': 0.70}
    else:
        cosmo_in = cosmo

    co_dist = cd.quick_distance_function(cd.comoving_distance, **cosmo_in)

    return co_dist(z)


def fit_mstar(mhalo):
    c = 0.129
    m0 = 10. ** 11.4
    alpha = 0.926
    beta = 0.261
    gama = 2.440

    out = c * mhalo * ((mhalo / m0) ** (-alpha) + (mhalo / m0)
                       ** beta) ** (-gama)
    return out


def halom_from_mstar(mstar):
    xbins = np.linspace(9.0, 16.0, 50)
    log_mstar = np.log10(mstar)
    y = np.log10(fit_mstar(10 ** xbins))
    func = interp1d(y, xbins)
    return func(log_mstar)


def matchidsorted(ids, targetid):
    """ Find id matches, return index in i1 that matches
    targetid; -1 if no match. """
    i1 = np.searchsorted(ids, targetid)
    if targetid == ids[i1]:
        ibest = i1
    else:
        ibest = -1
    return ibest


def acorr2(x, y, nbins=10, max_d=None, logbin=None):
    """ two points correlation function between x, y
        x is zip (d1, d2,..., dn), x is data array
        y is zip (d1, d2,..., dn), y is random points

        output:
            (bins, pairs counts)
        FIXME
    """
    pass


def matchids(id1, id2):
    """ Match two sets of ids.
        Returns:
          ibest -- array of indices of i1 that match i2; -1 if no match
    """
    id1 = np.atleast_1d(id1)
    id2 = np.atleast_1d(id2)
    indices = np.argsort(id1)
    idsorted = id1[indices]
    ibest = []
    id2_best = []
    for i in range(len(id2)):
        j = matchidsorted(idsorted, id2[i])
        if j >= 0:
            ibest += [indices[j]]
            id2_best += [i]

    return (np.array(ibest), np.array(id2_best))


def matchids_idx(id1, id2):
    """ Match two sets of ids.
        Returns:
          ibest -- array of indices of i1 that match i2; -1 if no match
    """
    id1 = np.atleast_1d(id1)
    id2 = np.atleast_1d(id2)
    indices = np.argsort(id1)
    idsorted = id1[indices]
    ibest = []
    id2_best = []
    for i in range(len(id2)):
        j = matchidsorted(idsorted, id2[i])
        if j >= 0:
            ibest += [indices[j]]
            id2_best += [i]

    return (np.array(ibest), np.array(id2_best))


def est_sm(g_minus_r, abs_r, redshift_band='0.0', band='r'):

    if band == 'r':
        mag_zero = 4.64

    gr = np.asarray(g_minus_r)
    abs_r = np.asanyarray(abs_r)
    lum = 10 ** (-0.4 * (abs_r - mag_zero))

    if redshift_band == '0.0':
        a = -0.7516
        b = 2.9838
        c = -0.5678

        sm_ratio = gr ** 2 * a + gr * b + c
        return sm_ratio * lum

    if redshift_band == '0.1':
        a = -0.6173
        b = 2.7839
        c = -0.42798

        sm_ratio = gr ** 2 * a + gr * b + c

        return sm_ratio * lum

def kde_hist_weight(data, xra, nbin=50, bandwidth=None, density=False,
                    weights=None, err=None, mirror=False, cdf=False):

    data = data[np.isfinite(data)]
    xmin, xmax = xra

    if mirror:
        idx = (data < xmin + 0.3)
        data = np.append(data, 2.0* xmin - data[idx])

    x_plot = np.linspace(xmin, xmax, nbin)
    kde_est = KDEUnivariate(data)
    fft_opt = False
    if weights is None:
        fft_opt = True
        weights_sum = len(data) * 1.0
    else:
        ftt_opt = False
        weights_sum = np.sum(weights)


    if bandwidth is not None:
        bw_in = bandwidth
    else:
        bw_in ='normal_reference'

    kde_est.fit(bw=bw_in, weights=weights, fft=fft_opt)
    if density:
        result = kde_est.evaluate(x_plot)
    else:
        result = kde_est.evaluate(x_plot) * weights_sum

    result_x = x_plot

    func = lambda x: kde_est.evaluate(x)

    if cdf:
        cdf = []
        for xx in x_plot:
            vv, _ = quad(func, xmin, xx)
            cdf.append(vv)

    if cdf:
        return result, result_x, np.array(cdf)
    else:
        return result, result_x


def kde_hist(data, xra, nbin=50, bandwidth=None, density=False,
             err=None, cdf=None, mirror=False):

    xmin, xmax = xra

    if mirror:
        # import ipdb; ipdb.set_trace()  # XXX BREAKPOINT

        idx = (data < xmin + 0.4)
        data = np.append(data, 2.0* xmin - data[idx])

        # idx = data > xmax - 0.3
        # data = np.append(data, 2.0*xmax - data[idx])

    data = data[:,np.newaxis]
    # import ipdb; ipdb.set_trace()  # XXX BREAKPOINT

    x_plot = np.linspace(xmin, xmax, nbin)[:, np.newaxis]
    kde = KernelDensity(bandwidth=bandwidth).fit(data)
    log_dens = kde.score_samples(x_plot)
    dens = np.exp(log_dens)
    if density:
        result_y = dens
    else:
        result_y = dens * len(data)
    result_x = x_plot

    func = lambda x: np.exp(kde.score(x))

    if cdf:
        cdf = []
        for xx in x_plot:
            vv, _ = quad(func, xmin, xx)
            cdf.append(vv)

    if err:
        nerr = 200
        boot_mean = np.zeros((nerr, nbin))
        for i in range(nerr):
            #s_data = kde.sample(n_samples=data.size)
            s_data = np.squeeze(bootstrap(data))
            offset = np.random.randn(s_data.size)*bandwidth
            s_data = s_data + offset

            kde_tmp = KDEUnivariate(s_data)
            kde_tmp.fit(bw=bandwidth )
            result = kde_tmp.evaluate(np.squeeze(x_plot)) * s_data.size

            #kde = KernelDensity(bandwidth=bandwidth).fit(s_data)
            #log_dens = kde.score_samples(x_plot) * s_data.size
            #result = np.exp(log_dens)
            boot_mean[i,:] = result

        boot_var = np.std(boot_mean, axis=0)
        return result_y, result_x , boot_var


    if cdf:
        return result_y, result_x, np.asarray(cdf)
    else:
        return result_y, result_x

if __name__ == "__main__":
    import doctest
    doctest.testmod()
