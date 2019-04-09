#! /usr/bin/env python2
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 qguo <qguo@lupus.aip.de>
#
# Distributed under terms of the MIT license.
#
# Last modified: 2015 Nov 05

"""
Modules for tools to find the nearest neighbour or
all neighbours in within certain diameter angle

"""
from __future__ import division
import numpy as np
from scipy.spatial import cKDTree as KDT
from scipy.spatial.distance import cdist as cdist


sph_old = False

class gals_tree:

    def __init__(self, ra_all, dec_all, degree=False):
        self.ra_all = ra_all
        self.dec_all = dec_all
        self.tree = build_sp_kdtree(self.ra_all, self.dec_all, degree=degree)


    def sph_query_ball_asp(self, ra_cen, dec_cen, r, degree=False):
        """
        Note: r: diamter angles seperation, not 3D distants
        """

        idx = sph_query_ball(self.tree, ra_cen, dec_cen, r, degree=degree)

        idx_result = []
        if np.array(r).size == 1:
            r_arr = np.zeros_like(ra_cen) + r
        else:
            r_arr = r

        ra_cen = np.atleast_1d(ra_cen)
        dec_cen = np.atleast_1d(dec_cen)
        r_arr = np.atleast_1d(r_arr)
        for idx_now, ra_now, dec_now, rnow in zip(idx, ra_cen, dec_cen,
                                            r_arr):
            angle_sp = sph_sp(self.ra_all[idx_now], self.dec_all[idx_now],
                                ra_now, dec_now, degree=degree)
            good_idx = np.squeeze(angle_sp) < rnow

            idx_this = np.array(idx_now)[np.atleast_1d(good_idx)]
            idx_result.append(idx_this)


        return idx_result


    def sph_query(self, ra_cen, dec_cen, degree=False):
        """
        Note: r: diamter angles seperation, not 3D distants
        """

        idx = self.tree.query(ra_cen, dec_cen, r, degree=degree)

        idx_result = []
        if np.array(r).size == 1:
            r_arr = np.zeros_like(ra_cen) + r
        else:
            r_arr = r

        ra_cen = np.atleast_1d(ra_cen)
        dec_cen = np.atleast_1d(dec_cen)
        r_arr = np.atleast_1d(r_arr)
        for idx_now, ra_now, dec_now, rnow in zip(idx, ra_cen, dec_cen,
                                            r_arr):
            angle_sp = sph_sp(self.ra_all[idx_now], self.dec_all[idx_now],
                                ra_now, dec_now, degree=degree)
            good_idx = np.squeeze(angle_sp) < rnow

            idx_this = np.array(idx_now)[np.atleast_1d(good_idx)]
            idx_result.append(idx_this)


        return idx_result


def radec2xyz(ra, dec, degree=False):

    ra = np.asarray(ra)
    dec = np.asarray(dec)

    if degree:
        ra_r = np.radians(ra)
        dec_r = np.radians(dec)
    else:
        ra_r = ra
        dec_r = dec

    x = np.cos(ra_r) * np.cos(dec_r)
    y = np.sin(ra_r) * np.cos(dec_r)
    z = np.sin(dec_r)
    coords = np.empty((x.size, 3))
    coords[:, 0] = x
    coords[:, 1] = y
    coords[:, 2] = z

    return x, y, z


def xyz2radec(x,y,z, degree=False):
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    vec = np.vstack((x,y,z)).T
    mod_vec = np.linalg.norm(vec, axis=1)

    dec = np.arcsin(z / mod_vec)
    ra = np.arccos(x / mod_vec / np.cos(dec))


    if x.ndim == 0:
        if y < 0:
            ra = np.pi * 2.0 - ra
    elif x.ndim > 0:
        idx = y < 0
        ra[idx] = np.pi * 2.0 - ra[idx]

    if degree:
        dec = dec /np.pi *180.
        ra = ra /np.pi *180.

    return ra,dec

def build_sp_kdtree(ra, dec, degree=False):
    """
    ra and dec are in radian, if not set degree to True
    """

    x, y, z = radec2xyz(ra, dec, degree=degree)
    coords = np.empty((x.size, 3))
    coords[:, 0] = x
    coords[:, 1] = y
    coords[:, 2] = z

    tree = KDT(coords, balanced_tree=False)

    return tree


def sph_query_ball_asp(ra_all, dec_all, ra_cen, dec_cen, r, degree=False):
    """
    Note: r: diamter angles seperation, not 3D distants
    """

    tree = build_sp_kdtree(ra_all, dec_all, degree=degree)

    idx = sph_query_ball(tree, ra_cen, dec_cen, r, degree=degree)
    # import ipdb; ipdb.set_trace()  # XXX BREAKPOINT



    idx_result = []
    if np.array(r).size == 1:
        r_arr = np.zeros_like(ra_cen) + r
    else:
        r_arr = r

    ra_cen = np.atleast_1d(ra_cen)
    dec_cen = np.atleast_1d(dec_cen)
    r_arr = np.atleast_1d(r_arr)
    for idx_now, ra_now, dec_now, rnow in zip(idx, ra_cen, dec_cen,
                                           r_arr):
        angle_sp = sph_sp(ra_all[idx_now], dec_all[idx_now],
                            ra_now, dec_now, degree=degree)

        good_idx = np.squeeze(angle_sp) < rnow

        # print np.squeeze(angle_sp)[good_idx]
        idx_this = np.array(idx_now)[np.atleast_1d(good_idx)]
        # print id_all[idx_now][np.atleast_1d(good_idx)]
        idx_result.append(idx_this)


    return idx_result, idx

def sph_query_ball_asp_with_tree(tree, ra_cen, dec_cen, r, degree=False):
    """
    Note: r: diamter angles seperation, not 3D distants
    """

    idx = sph_query_ball(tree, ra_cen, dec_cen, r, degree=degree)

    idx_result = []
    if np.array(r).size == 1:
        r_arr = np.zeros_like(ra_cen) + r
    else:
        r_arr = r

    for idx_now, ra_now, dec_now, rnow in zip(idx, ra_cen, dec_cen,
                                           r_arr):
        angle_sp = sph_sp(ra_all[idx_now], dec_all[idx_now],
                            ra_now, dec_now, degree=degree)
        good_idx = np.squeeze(angle_sp) < rnow

        idx_this = np.array(idx_now)[np.atleast_1d(good_idx)]
        idx_result.append(idx_this)


    return idx_result


def sph_query_ball(tree, ra, dec, r, degree=False):
    """
    Note: only query the the cubic distance, not diamter angles
          need to perform another selection
    """

    x, y, z = radec2xyz(ra, dec, degree=degree)

    coords = np.empty((x.size, 3))
    coords[:, 0] = x
    coords[:, 1] = y
    coords[:, 2] = z

    if np.array(r).size == 1:
        idx = tree.query_ball_point(coords, r)
    elif x.shape == np.array(r).shape:
        idx = []
        for pnow, rnow in zip(coords, r):
            idx_now = tree.query_ball_point(pnow, rnow)
            idx.append(idx_now)

    return idx

def sph_sp_xyz(x, y, z, ra2, dec2, degree=False):
    x1, y1, z1= x, y, z
    x2, y2, z2= radec2xyz(ra2, dec2, degree=degree)

    p1 = np.vstack((x1, y1, z1)).T
    p2 = np.vstack((x2, y2, z2)).T


    mod_p1 = np.linalg.norm(p1, axis=1)
    mod_p2 = np.linalg.norm(p2, axis=1)

    mod_p1p2= np.outer(mod_p1, mod_p2)


    dot = np.inner(p1, p2)/mod_p1p2
    print 'sph_sp_xyz may be wrong'
    import sys
    sys.exit()

    dot = np.clip(dot, -1.0, 1.0)

    angle = np.arccos(dot)

    return angle

def sph_sp(ra1, dec1, ra2, dec2, degree=False, return_diagnoal=False):
    """
    np.tensordot will create result (3x3) from input (3) and (3)
    Fixme: we may want from input [a1, a2, a3] and [b1, b2, b3]
    into [a1b1, a2b2, a3b3].

    So if ra1.size = 9,  ra2.size =1,  if create (9x1) result.
    So if ra1.size = 9,  ra2.size =9,  if create (9x9) result.
    But we may only want 9 results.
    """
    x1, y1, z1= radec2xyz(ra1, dec1, degree=degree)
    x2, y2, z2= radec2xyz(ra2, dec2, degree=degree)

    p1 = np.vstack((x1, y1, z1)).T
    p2 = np.vstack((x2, y2, z2)).T


    mod_p1 = np.linalg.norm(p1, axis=1)
    mod_p2 = np.linalg.norm(p2, axis=1)

    mod_p1p2= np.outer(mod_p1, mod_p2)


    if sph_old:
        dot = np.zeros(np.asarray(ra1).size)

        for ii, pp in enumerate(p1):
            dot[ii] = np.dot(pp, p2[0])/mod_p1p2[ii]
    else:

        dot = np.tensordot(p1, p2, axes=(1,1)) / mod_p1p2

    dot = np.clip(dot, -1.0, 1.0)

    angle = np.arccos(dot)

    if return_diagnoal:
        if dot.shape[0] == dot.shape[1]:
            return np.diagonal(angle)
        else:
            import sys
            print "cannot return diagonal, input a and a is not the same"
            sys.exit()
    else:
        return angle


def sph_sp_d(ra1, dec1, ra2, dec2, degree=False):
    #x1, y1, z1= radec2xyz(ra1, dec1, degree=degree)
    #x2, y2, z2= radec2xyz(ra2, dec2, degree=degree)

    ra1 = np.atleast_1d(ra1)
    dec1 = np.atleast_1d(dec1)

    if degree:
        ra1 = ra1 / 180. * np.pi
        dec1 = dec1 / 180. * np.pi

        ra2 = ra2 / 180. * np.pi
        dec2 = dec2 / 180. * np.pi

    cosangle = np.sin(dec1 + 0.5*np.pi) * np.sin(dec2 + 0.5*np.pi) * \
        np.cos(ra1 - ra2) + np.cos(dec1 + 0.5*np.pi) * \
        np.cos(dec2 + 0.5 *np.pi)

    cosangle = np.clip(cosangle, -1.0, 1.0)
    angle = np.arccos(cosangle)

    # print cosangle

    return angle

def cube_sp(ra1, dec1, ra2, dec2, degree=False):

    x1, y1, z1= radec2xyz(ra1, dec1, degree=degree)
    x2, y2, z2= radec2xyz(ra2, dec2, degree=degree)

    p1 = np.vstack((x1, y1, z1)).T
    p2 = np.vstack((x2, y2, z2)).T

    return cdist(p1, p2)

if __name__ == '__main__':


    ra = np.random.rand(1000) * np.pi
    dec = np.random.rand(1000) * np.pi - np.pi/2.0

    #ra0, dec0 = [0, 0.1, 0.2], [0, 0.1, 0.2]
    ra0, dec0 = [0.1, 0.1, 0.1], [-1.3, -1.3, -1.3]

    tree = build_sp_kdtree(ra, dec)

    idx = sph_query_ball(tree, ra0,  dec0, [1.8, 1.9, 2.0])

    for i, idx_now in enumerate(idx):

        print len(idx_now), 'kdetree: result, 3D distant'

        #print ra0, dec0
        #angle = sph_sp(ra[idx_now], dec[idx_now], ra0, dec0)
        #print angle

    idx_sp = sph_query_ball_asp(ra, dec, ra0, dec0, [1.8, 1.9, 2.0])


    for i, idx_now in enumerate(idx_sp):

        print len(idx_now), 'kdetree: result, angular SP'


    for ra_now, dec_now, r in zip(ra0, dec0, [1.8, 1.9, 2.0]):
        spnow = sph_sp(ra, dec, ra_now, dec_now)

        idx = spnow < r
        print ra[np.squeeze(idx)].shape, 'brutal resutls: angular SP'








