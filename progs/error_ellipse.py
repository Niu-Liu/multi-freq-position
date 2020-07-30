#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: error_ellipse.py
"""
Created on Mon Feb 24 10:45:59 2020

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from numpy import cos, sin, pi, deg2rad
from numpy.linalg import eig, inv
from scipy.stats import norm


__all__ = ["calc_ellipse_shape", "simulate_error_ellipse", "simulate_from_covariance"]


# -----------------------------  FUNCTIONS -----------------------------
def error_ellipse(ra_err, dec_err, ra_dec_corr, anticw=False):
    """Calculate the ovrall formal uncertainty.

    ovrall formal uncertainty = sqrt(RA_err^2+Dec_err^2+C*RA_err*Dec_err)

    Parameters
    ----------
    ra_err/dec_err : formal uncertainty of RA/Dec, usually in micro-as
    ra_dec_corr : correlation coeffient between RA and Dec, unitless.

    Returns
    ----------
    M : semi-major axis of the error ellipse
    m : semi-minor axis of the error ellipse
    pa : the position angle of the major axis of the error ellipse
    """

    cov = ra_dec_corr * ra_err * dec_err
    cov_mat = np.array([[ra_err**2, cov], [cov, dec_err**2]])

    from numpy.linalg import eig

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

    return M, m, pa


def rotation_z(x, y, pa, degree=False):
    """Rotate around Z-axis by a position angle.

    Parameters
    ----------
    x/y : ndarray of float
        original data
    pa : float
        positional angle, degree reckoned from x-axis and anti-clockwise.

    Returns
    -------
    x1/y1 : ndarray of float
        new data
    """

    if degree:
        alpha = deg2rad(pa)
    else:
        alpha = pa

    x1 = x * cos(alpha) - y * sin(alpha)
    y1 = x * sin(alpha) + y * cos(alpha)

    return x1, y1


def calc_ellipse_shape(M, m, pa, scale=1, start_from_xaxis=False):
    """Calculate the an ellipse shape curve.

    Parameters
    ----------
    M : float
        major axis
    m : floar
        minor axis
    pa : float
        positional angle, degree reckoned from x- or y-axis.

    Returns
    -------
    x1/y1 : ndarray of float
    """

    t = np.linspace(0, 2 * np.pi, 360)
    x = M * cos(t) * scale
    y = m * sin(t) * scale

    if start_from_xaxis:
        alpha = np.deg2rad(pa)
    else:
        alpha = np.deg2rad(90-pa)

    x1, y1 = rotation_z(x, y, alpha)

    return x1, y1


def ellipse_shape_from_cov(ra_err, dec_err, ra_dec_corr, scale=1,
                           start_from_xaxis=False):
    """Calculate the an ellipse shape curve.

    Parameters
    ----------
    M : float
        major axis
    m : floar
        minor axis
    pa : float
        positional angle, degree reckoned from x- or y-axis.

    Returns
    -------
    x1/y1 : ndarray of float
    """

    # Calculate Error ellipse parameters
    M, m, pa = error_ellipse(ra_err, dec_err, ra_dec_corr)

    # Simulate data points
    x, y = calc_ellipse_shape(M, m, pa, scale, start_from_xaxis)

    return x, y


def simulate_error_ellipse(M, m, pa, datasize=10000):
    """Simulate datapoints within an error ellipse.

    Parameters
    ----------
    M : float
        major axis
    m : floar
        minor axis
    pa : float
        positional angle, degree reckoned from x- or y-axis.

    Returns
    -------
    x/y : ndarray of float
        coordinate of datapoints
    """

    # Major axis
    norDist1 = norm(0, M)
    x = norDist1.rvs(datasize)

    # Minor axis
    norDist2 = norm(0, m)
    y = norDist2.rvs(datasize)

    x, y = rotation_z(x, y, 90-pa, degree=True)

    return x, y


def simulate_from_covariance(ra_err, dec_err, ra_dec_corr, datasize=10000):
    """Simulate datapoints within an error ellipse.

    Parameters
    ----------
    ra_err :
        formal uncertainty of RA(cos(Dec)), usually in micro-as
    dec_err :
        formal uncertainty of Dec, usually in micro-as
    ra_dec_corr :
        correlation coeffient between RA and Dec, unitless.

    Returns
    -------
    x/y : ndarray of float
        coordinate of datapoints
    """

    # Calculate Error ellipse parameters
    M, m, pa = error_ellipse(ra_err, dec_err, ra_dec_corr)

    # Simulate data points
    x, y = simulate_error_ellipse(M, m, pa, datasize)

    return x, y


def search_for_rho_peak(rho):
    """Find the peak offset
    """

    # A step of 30 deg
    rhomax = np.max(rho)
    bins = np.linspace(0, rhomax, 51)
    num, binedge = np.histogram(rho, bins)
    # Indice corresponding to maximum value
    ind = np.argmax(num)
    # get a rough estimate of the peak
    if 0 <= ind < 5:
        rhomin, rhomax = binedge[ind-1],  binedge[ind+1]
    elif ind == 0:
        rhomin, rhomax = -180,  binedge[1]
    else:
        rhomin, rhomax = binedge[ind],  180

    # Search in a small error
    bins = np.linspace(rhomin, rhomax, 11)
    num, binedge = np.histogram(rho, bins)
    # Indice corresponding to maximum value
    ind = np.argmax(num)
    rho_est = binedge[ind]

    # Divide the sample into two sample
    rho1 = np.sort(rho[rho < 0])
    rho2 = np.sort(rho[rho >= 0])
    nb1 = len(rho1)
    nb2 = len(rho2)

    # 1-sigma
    rate = (1-0.667) / 2
    nb = int(len(rho) * rate)

    if nb >= nb1:
        print("Could not estimate the down limit")
        rhodown = 0
    else:
        rhodown = rho1[nb1 - nb] + rho_est

    if nb >= nb2:
        print("Could not estimate the up limit")
        rhoup = rhomax
    else:
        rhoup = rho2[nb] + rho_est

    return rho_est, rhodown, rhoup


def pa_center(pa, pa0):
    """Wrap angles into the range of (pa0-180, pa0+180)
    """

    pa_diff = pa - pa0
    pa_diff = np.where(pa_diff < -180, pa_diff+360, pa_diff)
    pa_diff = np.where(pa_diff > 360, pa_diff-360, pa_diff)

    return pa_diff


def search_for_pa_peak(pa, pa0):
    """Find the peak angle
    """

    pa_cen0 = pa_center(pa, pa0)

    # A step of 30 deg
    bins = np.linspace(-180, 180, 51)
    num, binedge = np.histogram(pa_cen0, bins)
    # Indice corresponding to maximum value
    ind = np.argmax(num)
    # get a rough estimate of the peak
    if 0 <= ind < 5:
        pamin, pamax = binedge[ind-1],  binedge[ind+1]
    elif ind == 0:
        pamin, pamax = -180,  binedge[1]
    else:
        pamin, pamax = binedge[ind],  180

    # Search in a small error
    bins = np.linspace(pamin, pamax, 11)
    num, binedge = np.histogram(pa_cen0, bins)
    # Indice corresponding to maximum value
    ind = np.argmax(num)
    pa_est = binedge[ind] + pa0

    # Divide the sample into two sample
    pa_cen = pa_center(pa, pa_est)
    pa1 = np.sort(pa_cen[pa_cen < 0])
    pa2 = np.sort(pa_cen[pa_cen >= 0])
    nb1 = len(pa1)
    nb2 = len(pa2)

    # 1-sigma
    rate = (1-0.667) / 2
    nb = int(len(pa) * rate)

    if nb >= nb1:
        print("Could not estimate the down limit")
        padown = pa_est - 180
    else:
        padown = pa1[nb1 - nb] + pa_est

    if nb >= nb2:
        print("Could not estimate the up limit")
        paup = pa_est + 180
    else:
        paup = pa2[nb] + pa_est

    return pa_est, padown, paup
# --------------------------------- END --------------------------------
