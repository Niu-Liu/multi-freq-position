#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: calc_pa.py
"""
Created on Wed Nov  6 14:50:36 2019

@author: Neo(liuniu@smail.nju.edu.cn)

This program is written for calculating the position offset direction.

"""

import numpy as np


# -----------------------------  FUNCTIONS -----------------------------
def pa_calc(dra, ddec, anticw=False):
    """Calculate positional angle from positional offset.

    Parametes
    ---------
    dra : ndarray
        positional difference in R.A.(times cos(decl.))
    ddec : ndarray
        positional difference in declination

    Returns
    -------
    A : ndarray
        Angle (in degree) of positional offset vector towards to y-axis clockwisely
    """

    # Direction of positional offset vector from the y-axis count-clockwisely
    # in degree
    A = np.rad2deg(np.arctan2(dra, ddec))  # clockwise

    if anticw:
        # anticlockwise
        A = np.where(A < 0, -A, 360 - A)
    else:
        # clockwise
        A = np.where(A < 0, 360 + A, A)

    return A


def inclu_angle_calc(dra1, ddec1, dra2, ddec2):
    """Calculate included angle between two positional offset vectors.


    Parametes
    ---------
    dra1, ddec1: ndarray
        components of the fisrt positional offset vector
    dra2, ddec2: ndarray
        components of the second positional offset vectors

    Returns
    -------
    ang : ndarray
        Angle (in degree) between positional offset vectors
    """

    cosine = dra1 * dra2 + ddec1 * ddec2
    mod1 = np.sqrt(dra1**2 + ddec1**2)
    mod2 = np.sqrt(dra2**2 + ddec2**2)

    ang = np.rad2deg(np.arccos(cosine/mod1/mod2))

    return ang


def eepa_calc_single(ra_err, dec_err, ra_dec_corr, anticw=False):
    """Calculate the ovrall formal uncertainty.

    ovrall formal uncertainty = sqrt(RA_err^2+Dec_err^2+C*RA_err*Dec_err)

    Parameters
    ----------
    ra_err/dec_err : formal uncertainty of RA/Dec, usually in micro-as
    ra_dec_corr : correlation coeffient between RA and Dec, unitless.

    Returns
    ----------
    pa : the position angle of the major axis of the error ellipse
    """

    from numpy.linalg import eig

    cov = ra_dec_corr * ra_err * dec_err
    cov_mat = np.array([[ra_err**2, cov], [cov, dec_err**2]])

    eig_val, eig_vec = eig(cov_mat)

    if eig_val[0] > eig_val[1]:
        M2 = eig_val[0]
        m2 = eig_val[1]
        vec_M = eig_vec[:, 0]
    else:
        M2 = eig_val[1]
        m2 = eig_val[0]
        vec_M = eig_vec[:, 1]

    M, m = np.sqrt(M2), np.sqrt(m2)

    # Calculate the positi angle counted anti-colockwise from the declination axis
    if vec_M[1] == 0:
        # It rarely happens.
        if vec_M[0] > 0:
            pa0 = 90
        else:
            pa0 = -90
    else:
        pa0 = np.rad2deg(np.arctan(vec_M[0]/vec_M[1]))

    if anticw:
        if pa0 <= 0:
            pa = -pa0
        else:
            pa = 180 - pa0
    else:
        if pa0 <= 0:
            pa = 180 + pa0
        else:
            pa = pa0

    return pa


def eepa_calc(ra_err, dec_err, ra_dec_corr, anticw=False):
    """Calculate the ovrall formal uncertainty.

    ovrall formal uncertainty = sqrt(RA_err^2+Dec_err^2+C*RA_err*Dec_err)

    Parameters
    ----------
    ra_err/dec_err : formal uncertainty of RA/Dec, usually in micro-as
    ra_dec_corr : correlation coeffient between RA and Dec, unitless.

    Returns
    ----------
    pa : the position angle of the major axis of the error ellipse
    """

    pa = np.zeros_like(ra_err)
    for i in range(len(ra_err)):
        pa[i] = eepa_calc_single(ra_err[i], dec_err[i], ra_dec_corr[i], anticw)

    return pa


def wrap_angle(pa):
    """Wrap angle > 180 over 0-180

    Parameter
    ---------
    pa : ndarray
        position angle in unit of degree over 0-360

    Return
    ------
    pa1 : ndarray
        warpped position angle in unit of degree over 0-360
    """

    pa = np.array(pa)

    pa = np.where(pa > 180, pa - 360, pa)
    pa = np.where(pa < -180, 360 + pa, pa)

    return pa
# --------------------------------- END --------------------------------
