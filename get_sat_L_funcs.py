#! /usr/bin/env python2
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 dccf87 <dccf87@gqlinux>
#
# Distributed under terms of the MIT license.
#
# Last modified: 2019 Apr 08

"""
cacluate the satellite luminosity funcitons
and number density profiles.
"""


import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from gqcos import lum_distf, abs_mf, kde_hist, kde_hist_weight, bootstrap
from spdist import sph_query_ball_asp, sph_sp, gals_tree
from gqplot import logerr
from data_input import get_sdss_data
import sys

def sat_weigh(mj, mlimit, binw):
    mj = np.atleast_1d(mj)
    mlimit = np.atleast_1d(mlimit)

    w = np.zeros_like(mj)

    idx = mj + binw < mlimit
    if np.any(idx):
        w[idx] = 1.0

    idx = (mj + binw > mlimit) & (mj < mlimit)
    if np.any(idx):
        w[idx] = 1.0

    return w

def prim_weigh(mj, mlimit, binw):
    mj = np.atleast_1d(mj)
    mlimit = np.asarray(mlimit)

    w = np.zeros_like(mj)

    idx = mj + binw < mlimit
    if np.any(idx):
        w[idx] = 1.0

    idx = (mj + binw > mlimit) & (mj < mlimit)
    if np.any(idx):
        w[idx] = (mlimit - mj[idx]) / binw

    return w


job_id_all = ['test']
mc_all = [-21., -22., -23.]
# mc_all = [-21.,-23.]
color_all = ['k', 'b', 'r']
# mc_all = [-21.]

data_s_all = ['SDSS']
prim_files = ['prims.h5']
neigs_files = ['neigs.h5']
xmin,xmax = [0.0, 10.0]
nbin = 20
old_method = False
kde_method = True
bootstrap_err = True



for data_s, p_file, n_file in zip(data_s_all, prim_files, neigs_files):

    of = h5.File('../output/' + data_s + '_plot.h5', 'w')

    for job_id, data_s in zip(job_id_all, data_s_all):
        fin = '../output/' + data_s + '_' + job_id + '.h5'
        # print 'unpacking data started'
        # with gzip.GzipFile(fin,'r') as fsat:
            # data = ujson.load(fsat)
        # print 'unpacking data fnished'

        data = h5.File(fin, 'r')


        if 'SDSS' in data_s:
            fp = h5.File('../catas/' + p_file, mode='r')
            fn = h5.File('../catas/' + n_file, mode='r')

            cosmo = {}
            for item in data['cosmo']:
                cosmo[item] = np.float(data['cosmo'][item][...])

            print 'Cosmology', cosmo

            band = str(data['-21.0']['band'][...])
            kcorr = str(data['-21.0']['kcorr'][...])
            pgals, ngals = get_sdss_data(fp, fn, band, kcorr, **cosmo)
            print 'read data files'


        for ii, mc in enumerate(mc_all):
            mc = str(mc)
            print 'mc: ', mc
            prim_id = data[mc]['prim_id'][...]
            prim_gals = pgals[prim_id]
            inner_id = data[mc]['inner_id'][...]
            outer_id = data[mc]['outer_id'][...]
            mag_limit = data[mc]['mag_limit'][...]
            mask_i = data[mc]['mask_i'][...]
            mask_o = data[mc]['mask_o'][...]


            r_i = data[mc]['r_i'][...]
            r_o = data[mc]['r_o'][...]
            area_f = r_i ** 2 / (r_o ** 2 - r_i ** 2)


            if old_method:
                sat_fg_dm = []
                sat_bg_dm = []

                for p, fg_id, bg_id in zip(prim_gals, inner_id, outer_id):

                    fgs = ngals[fg_id]
                    bgs = ngals[bg_id]

                    fg_dm = fgs['mag'] - p['band_mag']
                    bg_dm = bgs['mag'] - p['band_mag']
                    sat_fg_dm.extend(fg_dm)
                    sat_bg_dm.extend(bg_dm)
                        # import ipdb; ipdb.set_trace()  # XXX BREAKPOINT


                hy_fg, hx = np.histogram(sat_fg_dm, nbin, range=[xmin, xmax])
                hy_bg, hx = np.histogram(sat_bg_dm, nbin, range=[xmin, xmax])

                binw = hx[1] - hx[0]

                px = []
                py = []
                num_p = []

                for i in range(hx[:-1].size):
                    # idx = (mag_limit - prim_gals['band_mag'] > hx[i])
                    idx = (mag_limit - prim_gals['band_mag'] > hx[i] + 0.5*binw)
                    # idx = (prim_gals['mag'] + hx[i] + binw < mag_limit)
                    if np.any(idx):
                        # prim_cnt = np.sum(prim_weigh(prim_gals['mag'] + hx[i],
                                        # mag_limit, binw))

                        prim_cnt = np.float(prim_gals['band_mag'][idx].size)


                        sat_num = (hy_fg[i] - hy_bg[i] * area_f)/prim_cnt/binw
                        py.append(sat_num)
                        px.append(hx[i] + 0.5*binw)
                        num_p.append(prim_cnt)
                        print hx[i] + 0.5*binw, sat_num, prim_cnt, prim_cnt

                        if i == 25:
                            import ipdb; ipdb.set_trace()  # XXX BREAKPOINT


                plt.plot(px, np.log10(py), 'k--', color =color_all[ii], lw=2.0, ms=5.0)

            if kde_method:

                sat_fg_dm = []
                sat_bg_dm = []
                delta_m = mag_limit - prim_gals['band_mag']

                for p, fg_id, bg_id, m_i, m_o in zip(prim_gals, inner_id,
                                                     outer_id, mask_i, mask_o):

                    fgs = ngals[fg_id]
                    bgs = ngals[bg_id]

                    fg_dm = fgs['mag'] - p['band_mag']
                    bg_dm = bgs['mag'] - p['band_mag']
                    # bg_weights = np.zeros_like(bg_dm)
                    # fg_weights = np.zeros_like(fg_dm)
                    # fg_weights += 1.0/m_i
                    # bg_weights += 1.0/m_o
                    # sat_fg_weights = np.append(sat_fg_weights, fg_weights)
                    # sat_bg_weights = np.append(sat_bg_weights, bg_weights)
                    sat_fg_dm.extend(fg_dm)
                    sat_bg_dm.extend(bg_dm)



                dmax = min(np.array(sat_fg_dm).max(), np.array(sat_fg_dm).max())-0.3
                print dmax

                if mc == '-21.0':
                    bw = 0.4
                else:
                    bw = 0.3
                sat_fg_py, px = kde_hist_weight(np.asarray(sat_fg_dm), [0.5, dmax],
                                         bandwidth=bw, mirror=True)
                sat_bg_py, px = kde_hist_weight(np.asarray(sat_bg_dm), [0.5, dmax],
                                         bandwidth=bw, mirror=True)
                _, _, cdf= kde_hist_weight(delta_m,[0.5, dmax], bandwidth=bw, cdf=True)

                # plt.plot(px, np.log10((1-cdf)*delta_m.size), 'k-', lw=2.0)
                # plt.plot(px, np.log10((sat_fg_py- area_f*sat_bg_py)/((1-cdf)*delta_m.size)), 'k-', lw=2.0)


                # plt.show()
                # import sys
                # sys.exit()

                lf_px = px
                lf_py = (sat_fg_py- area_f*sat_bg_py)/((1-cdf)*delta_m.size)

                if bootstrap_err:
                    nerr = 5
                    boot_mean = np.zeros((nerr, px.size))

                    for i in range(nerr):

                        sat_fg_dm = []
                        sat_bg_dm = []

                        print 'bootstrap, time i'
                        boot_idx = bootstrap(prim_gals, return_index=True)

                        s_prim_gals = prim_gals[boot_idx]
                        s_inner_id = inner_id[boot_idx]
                        s_outer_id = outer_id[boot_idx]
                        delta_m = mag_limit - s_prim_gals['band_mag']

                        for p, fg_id, bg_id in zip(prim_gals, inner_id, outer_id,):

                            fgs = ngals[fg_id]
                            bgs = ngals[bg_id]

                            fg_dm = fgs['mag'] - p['band_mag']
                            bg_dm = bgs['mag'] - p['band_mag']
                            sat_fg_dm.extend(fg_dm)
                            sat_bg_dm.extend(bg_dm)


                        sat_fg_dm = np.array(sat_fg_dm)
                        sat_bg_dm = np.array(sat_bg_dm)

                        dmax = min(np.array(sat_fg_dm).max(), np.array(sat_fg_dm).max())-0.3

                        bandwidth = 0.3

                        offset = np.random.randn(sat_fg_dm.size)*bandwidth
                        sat_fg_dm += offset

                        sat_fg_py, px = kde_hist_weight(np.asarray(sat_fg_dm), [0.5, dmax],
                                                bandwidth=0.3, mirror=True)

                        offset = np.random.randn(sat_bg_dm.size)*bandwidth
                        sat_bg_dm += offset

                        sat_bg_py, px = kde_hist_weight(np.asarray(sat_bg_dm), [0.5, dmax],
                                                bandwidth=0.3, mirror=True)

                        bandwidth = 0.2
                        offset = np.random.randn(delta_m.size)*bandwidth
                        delta_m += offset
                        py, px, cdf= kde_hist_weight(delta_m,[0.5, dmax], bandwidth=0.2, cdf=True)

                        boot_lf_py = (sat_fg_py- area_f*sat_bg_py)/((1-cdf)*delta_m.size)

                        #kde = KernelDensity(bandwidth=bandwidth).fit(s_data)
                        #log_dens = kde.score_samples(x_plot) * s_data.size
                        #result = np.exp(log_dens)
                        boot_mean[i,:] = boot_lf_py

                    boot_err = np.std(boot_mean, axis=0)

                    output = np.percentile(boot_mean,[5, 95], axis=0)
                    # import ipdb; ipdb.set_trace()  # XXX BREAKPOINT

                    errors = logerr(lf_py, boot_err)

                    plt.fill_between(px, np.log10(lf_py) - errors[0],
                                    np.log10(lf_py)+errors[1], color=color_all[ii],
                                    alpha=0.5, edgecolor='none')

                plt.plot(lf_px, np.log10(lf_py), '-', color=color_all[ii],
                         lw=2.0, alpha=0.8, label='$M_c = $' + mc)


                grp = of.create_group(str(mc))
                grp.create_dataset('lf_px', data=lf_px)
                grp.create_dataset('lf_py', data=lf_py)
                if bootstrap_err:
                    grp.create_dataset('errs', data=boot_err)
                    grp.create_dataset('p5', data=output[0])
                    grp.create_dataset('p95', data=output[1])

                # method 3



            if mc == '-21.0':
                px_o,py_o,err = np.loadtxt('../../sat/results/21.00dr8_band_r_clean_nei_hehe_dm_bin_0.50_dm_f_0.50_ap_2.50_dzs_600.00_ri2_0.30_ri1_0.00_ro2_0.60_ro1_0.30_magl_20.50_final_sat.dat_lf_all10.00',
                                           unpack=True, usecols=(0,1,2))
                plt.plot(px_o, np.log10(py_o), 'ko', ms=5.0, alpha=0.8)


            if mc == '-22.0':
                px_o,py_o,err = np.loadtxt('../../sat/results/22.00dr8_band_r_clean_nei_dm_bin_0.50_dm_f_0.50_ap_2.50_dzs_1000.00_ri2_0.40_ri1_0.00_ro2_0.80_ro1_0.40_magl_20.50_final_sat.dat_lf_all10.00',
                                           unpack=True, usecols=(0,1,2))
                plt.plot(px_o, np.log10(py_o), 'bo', ms=5.0, alpha=0.8,
                         label='histogram')

            if mc == '-23.0':
                px_o,py_o,err = np.loadtxt('../../sat/results/23.00dr8_band_r_clean_nei_hehe_dm_bin_0.50_dm_f_0.50_ap_2.50_dzs_1200.00_ri2_0.55_ri1_0.00_ro2_0.90_ro1_0.55_magl_20.50_final_sat.dat_lf_all10.00',
                                           unpack=True, usecols=(0,1,2))
                plt.plot(px_o, np.log10(py_o), 'ro', ms=5.0, alpha=0.8)


        plt.xlabel(r'$\Delta M$')
        plt.ylabel(r'$\log(dN/dM)$')
        plt.ylim([-2, 1.5])
        plt.legend(frameon=False, loc='best')
        plt.tight_layout()
        plt.show()
        # plt.savefig('test.pdf')


        of.close()
        fp.close()
        fn.close()
