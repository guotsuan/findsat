#! /usr/bin/env python2
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 dccf87 <dccf87@gqlinux>
#
# Distributed under terms of the MIT license.
#
# Last modified: 2017 Jan 30

"""
get different data set
"""
import numpy as np
from gqcos import lum_distf, abs_mf
from spdist import sph_sp
import sys
import mask_tools


def get_sdss_data(fp, fn,  band='r', kcorr='00', **cosmo):

    pgal_type = [('ra', np.float), ('dec', np.float),
                 ('redshift', np.float), ('objid', np.int64),
                 ('dist', np.float), ('sn', np.int),
                 ('mag', np.float), ('band_mag', np.float),
                 ('kr', np.float), ('abs_m', np.float),
                 ('color', np.float), ('sm', np.float)]


    ngal_type = [('ra', np.float), ('dec', np.float),
                 ('flags', np.int64), ('objid', np.int64),
                 ('mag', np.float), ('spz', np.float),
                 ('phoz', np.float), ('zerr', np.float),
                 ('sn', np.int), ('mag_g', np.float)]


    # primary galaxies
    p_ra = fp['RA'][...]
    p_dec = fp['DEC'][...]
    p_redshift = fp['RS'][...]
    p_redshift = np.clip(p_redshift, 1e-7, 5.0)
    p_dist = lum_distf(p_redshift, **cosmo) / (1. + p_redshift) ** 2
    p_objid = fp['OBJID'][...]
    p_color = fp['COLOR'][...]
    p_sm = fp['K_SM_00'][...]
    p_seq = np.arange(p_ra.size)


    if '00' in kcorr:
        kcorr_red = '_00_'
    elif '01' in kcorr:
        kcorr_red = '_01_'
    else:
        print 'kcorr setting error, only 00 or 01 for now'
        sys.exit(1)

    if 'r' in band:
        band_name ='R'
        kcorr_name = 'KCORR' + kcorr_red + band_name
        print 'kcorr_name', kcorr_name

        kcorr_raw = fp[kcorr_name][...]
        good = np.isfinite(kcorr_raw)
        kcorr_raw[~good] = 9999.
        p_mag = fp[band_name][...] - kcorr_raw
        p_band = fp[band_name][...]
    else:
        print 'band selecting error:'
        sys.exit(1)

    p_abs_m = abs_mf(p_mag, p_redshift, **cosmo)
    prim_gals = np.zeros_like(p_ra, pgal_type)

    prim_gals['ra'] = p_ra
    prim_gals['dec'] = p_dec
    prim_gals['redshift'] = p_redshift
    prim_gals['dist'] = p_dist
    prim_gals['objid'] = p_objid
    prim_gals['sn'] = p_seq
    prim_gals['kr'] = kcorr_raw
    prim_gals['mag'] = p_mag
    prim_gals['abs_m'] = p_abs_m
    prim_gals['band_mag'] = p_band
    prim_gals['color'] = p_color
    prim_gals['sm'] = p_sm


    # neighboring galaxies
    n_ra = fn['RA'][...]
    n_dec = fn['DEC'][...]
    n_flags = fn['FLAGS'][...]
    n_objid = fn['OBJID'][...]
    n_mag = fn[band_name][...]
    n_spz = fn['SPZ'][...]
    n_phoz = fn['PHOZ'][...]
    n_zerr= fn['ZERR'][...]
    n_seq = np.arange(n_ra.size)
    n_mag_g = fn['G'][...]

    ngals = np.zeros_like(n_ra, ngal_type)
    ngals['ra'] = n_ra
    ngals['dec'] = n_dec
    ngals['flags'] = n_flags
    ngals['mag'] = n_mag
    ngals['spz'] = n_spz
    ngals['phoz'] = n_phoz
    ngals['zerr'] = n_zerr
    ngals['sn'] = n_seq
    ngals['objid'] = n_objid
    ngals['mag_g'] = n_mag_g

    return prim_gals, ngals


def sdss_mask(pra, pdec, r_i, r_o):
    """
    par, pdec, are the ra, dec for querying
    r_i, r_o are the angular sepration radius
    """

    num = 100

    radius = 1.5 * r_o
    ran_ra = (radius  - 2.0 * radius * np.random.rand(3500)) / np.cos(pdec)
    ran_dec = (radius - 2.0 * radius * np.random.rand(3500))

    ra_list = pra + ran_ra
    dec_list = pdec + ran_dec

    sp = sph_sp(ra_list, dec_list, pra, pdec)

    idx_i = np.squeeze(sp) < r_i

    idx_o = (np.squeeze(sp) > r_i) & (np.squeeze(sp) < r_o )

    ran_ra_s = np.random.choice(ra_list[idx_i], min(num, ra_list[idx_i].size), replace=False)
    ran_dec_s = np.random.choice(dec_list[idx_i], min(num, dec_list[idx_i].size), replace=False)

    result_i = 0.0

    for ra, dec in zip(ran_ra_s, ran_dec_s):
        compl = mask_tools.query_mask(ra, dec)
        if compl > 0:
            compl = 1.0
        else:
            compl = 0.0

        result_i += compl

    result_i /= np.float(num)

    ran_ra_s = np.random.choice(ra_list[idx_o], min(num, ra_list[idx_o].size), replace=False)
    ran_dec_s = np.random.choice(dec_list[idx_o], min(num, dec_list[idx_o].size), replace=False)

    result_o = 0.0
    for ra, dec in zip(ran_ra_s, ran_dec_s):
        compl = mask_tools.query_mask(ra, dec)
        if compl > 0:
            compl = 1.0
        else:
            compl = 0.0


        result_o += compl

    result_o /= num
    return (result_i, result_o)




