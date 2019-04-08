#! /usr/bin/env python2
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 dccf87 <dccf87@gqlinux>
#
# Distributed under terms of the MIT license.
#
# Last modified: 2015 Nov 11

"""
find the isolated priamries in redshifts space
"""

import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from gqcos import lum_distf, abs_mf
from spdist import sph_query_ball_asp, sph_sp, gals_tree
from data_input import get_sdss_data,sdss_mask
import sys
import ujson
import gzip
import cPickle as pickle


outdir = '../output/'
meta_info = {}
# all units are without h0 by default

output_mask = True

h0 = 0.70
factor = -5.0 * np.log10(h0)
kcorr = '00'
band = 'r'
dm_f = 0.5
iso_dzs = 600.
min_phoz_err=0.05
a_p = 2.5
mag_limit = 20.5
mask_limit = 0.7
vel_c = 3.e5 # km/s

bad_flag = 0x40000 | 0x80000000000

#parameters
# bin width for prim selection
dm_bin = 0.5
cosmo = {'omega_M_0': 0.3, 'omega_lambda_0': 0.7, 'omega_k_0': 0.0, 'h':h0 }
data_s = 'SDSS'

M_c_all = [-21.,-22., -23.]
# M_c_all = [-21.]
R_i_all = [0.3, 0.4, 0.55]
R_o_all = [0.6, 0.8, 0.9]
R_iso_all = [0.6, 0.8, 1.1]
dzs_all = [600., 1000., 1200]



prim_files = ['prims.h5']
neigs_files = ['neigs.h5']

id_for_jobs =['test']


for p_file, n_file, job_id in zip(prim_files, neigs_files, id_for_jobs):
    print '../catas/' + p_file
    fp = h5.File('../catas/' + p_file, mode='r')
    fn = h5.File('../catas/' + n_file, mode='r')


    if 'SDSS' in data_s:
        pgals, ngals = get_sdss_data(fp, fn, band = 'r', kcorr = '00', **cosmo)

        # only use galaxies brighter than mag_limit
        ngals = ngals[ngals['mag'] < mag_limit]

        # only use primary galaxies in north cap


    if max(pgals['ra']) > 2 * np.pi:
        print 'I guest RA is in dgree, converting'
        pgals['ra'] = np.radians(pgals['ra'])
        pgals['dec'] = np.radians(pgals['dec'])

        ngals['ra'] = np.radians(ngals['ra'])
        ngals['dec'] = np.radians(ngals['dec'])
    else:
        print 'RA and DEC seems in unit of radians, continue...'

    print 'buding trees..'
    gals_tree = gals_tree(ngals['ra'], ngals['dec'])
    print 'tree is  built'

    for nn, mc in enumerate(M_c_all):
        # only north cap
        print 'start to search: ', mc
        p_idx = (pgals['abs_m'] > mc - dm_bin) & \
            (pgals['abs_m']< mc + dm_bin) & \
            (pgals['ra'] > np.radians(110.)) & \
            (pgals['ra'] < np.radians(290.)) & \
            (pgals['band_mag'] < 18.5)

        a_sp = max(R_iso_all[nn], R_o_all[nn]) / pgals['dist'][p_idx]
        close_idx = gals_tree.sph_query_ball_asp(pgals['ra'][p_idx],
                                       pgals['dec'][p_idx], a_sp)

        print 'find neighbours finished'

        #p_idx = (pgals['objid'] == 1237663788487278839 ) | \
              #(pgals['objid'] == 1237663917872054341 )

        #p_idx = (pgals['objid'] == 1237663655343161873 ) | \
              #(pgals['objid'] == 1237673706650534285)

        #p_idx = (pgals['objid'] == 1237663655880294818) | \
              #(pgals['objid'] == 1237673704503116115)

        #p_idx = (pgals['objid'] == 1237659163350663684)
        #p_idx = (p_objid == 1237665026517958779)

        #p_idx = (p_objid == 1237663789023560204) | \
            #(p_objid == 1237663789023560233) | \
            #(p_objid == 1237663789023756776) | \
            #(p_objid == 1237663789023560111)


        #p_idx = (p_abs_m > mc - dm_bin) & (p_abs_m < mc + dm_bin) \
            #& (p_ra > np.radians(110.)) & (p_ra < np.radians(290.)) \
            #& (p_mag -  < 18.5) & (p_redshift > 0.0 )

        print 'primary candiats', mc, pgals[p_idx].size

        p_good_sn = pgals['sn'][p_idx]

        good_mask = np.ones_like(p_good_sn, np.bool)
        p_mask_i= np.zeros_like(p_good_sn, np.float)
        p_mask_o= np.zeros_like(p_good_sn, np.float)

        giter = np.nditer([p_good_sn, good_mask, p_mask_i, p_mask_o],
                          op_flags=[['readonly'], ['readwrite'],
                                    ['readwrite'], ['readwrite']])

        inner_ids = []
        outer_ids = []

        for (p, mask_now, mask_i, mask_o), idx_this in zip(giter, close_idx):

            if output_mask:
                mask_i[...], mask_o[...] = sdss_mask(pgals['ra'][p], pgals['dec'][p],
                                           R_i_all[nn]/pgals['dist'][p],
                                           R_o_all[nn]/pgals['dist'][p])

                if mask_i < mask_limit or mask_o < mask_limit:
                    mask_now[...] = False
            else:
                # no cut in completness, every primary is included.
                mask_now[...] = True

            if len(idx_this) == 1 :
                if pgals['objid'][p] == ngals['objid'][idx_this[0]]:
                    pass
                else:
                    print 'theretical impossible, check code'
                    sys.exit()

            else:
                sp = np.squeeze(sph_sp(ngals['ra'][idx_this], ngals['dec'][idx_this],
                        pgals['ra'][p], pgals['dec'][p]))

                # find all bring neighbours which is not prim itself.
                bright_idx = (ngals['mag'][idx_this] -
                              pgals['band_mag'][p] <= dm_f) & \
                    (ngals['objid'][idx_this] != pgals['objid'][p])

                iso_close = (sp * pgals['dist'][p] < R_iso_all[nn])


                if np.any(bright_idx):
                    spz_idx = ngals['spz'][idx_this][bright_idx & iso_close] > -100
                    bright_spz = ngals['spz'][idx_this][bright_idx & iso_close][spz_idx]

                    if np.any(spz_idx):
                        vdiff = np.abs(bright_spz - pgals['redshift'][p]) * vel_c
                        vdiff = vdiff / (1+ 0.5*(bright_spz + pgals['redshift'][p]))
                        if np.any(vdiff < iso_dzs):
                            mask_now[...] = False

                    if mask_now[...]:
                        phoz_idx = ~spz_idx
                        bright_phoz = ngals['phoz']\
                            [idx_this][bright_idx & iso_close][phoz_idx]
                        bright_phoz_err = ngals['zerr'] \
                            [idx_this][bright_idx & iso_close][phoz_idx]

                        if np.any(phoz_idx):
                            phoerr_th = a_p * np.clip(bright_phoz_err,
                                                      min_phoz_err, 1.e5)
                            if np.any(np.abs(bright_phoz -
                                             pgals['redshift'][p]) < phoerr_th):
                                mask_now[...] = False


            # let us select inner and outer galaxies

            if mask_now:
                close = (ngals['mag'][idx_this] - pgals['band_mag'][p] > dm_f) & \
                    (ngals['objid'][idx_this] != pgals['objid'][p]) & \
                    (ngals['mag'][idx_this] < mag_limit) & \
                    (ngals['flags'][idx_this] & bad_flag == 0)


                # if have spz

                vdiff = np.abs(ngals['spz'][idx_this] - pgals['redshift'][p]) * vel_c
                vdiff = vdiff / (1+ 0.5*(ngals['spz'][idx_this] + pgals['redshift'][p]))
                spz_close = (ngals['spz'][idx_this] > -100 ) & \
                            (vdiff < dzs_all[nn] )

                # if no spz, let ask phoz for help

                phoz_err = ngals['zerr'][idx_this]
                phoerr_th = a_p * np.clip(phoz_err, min_phoz_err, 1.e5)
                phoz_close = (ngals['spz'][idx_this] < -100 ) & \
                    (np.abs(ngals['phoz'][idx_this] - pgals['redshift'][p]) < phoerr_th)


                inner = (sp * pgals['dist'][p] < R_i_all[nn])
                outer = (sp * pgals['dist'][p] >= R_i_all[nn]) & \
                        (sp * pgals['dist'][p] < R_o_all[nn])

                inner_gals = inner & close & (spz_close | phoz_close)
                outer_gals = outer & close & (spz_close | phoz_close)



                inner_ids.append(ngals['sn'][idx_this][inner_gals])
                outer_ids.append(ngals['sn'][idx_this][outer_gals])
            #print 'inner number'
            #print ngals['mag'][idx_this][inner_gals].size
            #print ngals['mag'][idx_this][outer_gals].size

        meta_info[mc] = {}
        # meta_info[mc]['cosmo'] = cosmo
        meta_info[mc]['dm_f'] = dm_f
        meta_info[mc]['dm_bin'] = dm_f
        meta_info[mc]['kcorr'] = kcorr
        meta_info[mc]['band'] = band
        meta_info[mc]['iso_dzs'] = iso_dzs
        meta_info[mc]['min_phoz_err'] = min_phoz_err
        meta_info[mc]['a_p'] = a_p
        meta_info[mc]['mag_limit'] = mag_limit
        meta_info[mc]['prim_id'] = list(p_good_sn[good_mask])
        meta_info[mc]['inner_id'] = inner_ids
        meta_info[mc]['outer_id'] = outer_ids
        meta_info[mc]['r_i'] = R_i_all[nn]
        meta_info[mc]['r_o'] = R_o_all[nn]
        meta_info[mc]['dzs'] = dzs_all[nn]
        meta_info[mc]['data_s'] = data_s
        meta_info[mc]['mask_limit'] = mask_limit
        meta_info[mc]['mask_i'] = p_mask_i[good_mask]
        meta_info[mc]['mask_o'] = p_mask_o[good_mask]



        print mc, 'finished', 'valid prims: ', p_good_sn[good_mask].size, \
               len(inner_ids)

    # save data
    fsat = data_s + '_' + job_id + '.h5'
    print 'outfile: ', fsat

    #output cosmology
    fo = h5.File(outdir + fsat, 'w')
    g_cosmo = fo.create_group('cosmo')
    for key in cosmo:
        g_cosmo.create_dataset(key, data=cosmo[key])
        print key, cosmo[key]


    for mc in M_c_all:
        grp = fo.create_group(str(mc))
        print 'group: ', str(mc)

        for key in meta_info[mc]:
            if key == 'inner_id' or key == 'outer_id':
                # print key
                dt = h5.special_dtype(vlen=np.int)
                # import ipdb; ipdb.set_trace()  # XXX BREAKPOINT

                grp.create_dataset(key, data=np.asarray(meta_info[mc][key]), dtype=dt)
            elif key == 'prim_id':
                grp.create_dataset(key, data=np.squeeze(meta_info[mc][key]))
            else:
                grp.create_dataset(key, data=meta_info[mc][key])



    fo.close()
    fp.close()
    fn.close()
