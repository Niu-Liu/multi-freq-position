#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: pos_simu.py
"""
Created on Mon Feb 24 13:31:14 2020

@author: Neo(liuniu@smail.nju.edu.cn)

Simulate the data points considering the covariance.

"""

from astropy.table import Table, vstack
import numpy as np
import time

from error_ellipse import simulate_from_covariance
from error_ellipse import search_for_pa_peak, search_for_rho_peak
from calc_pa import pa_calc, wrap_angle


# -----------------------------  FUNCTIONS ---------------------------------
def calc_sim_offset(comsou, i):
    """
    """
    souname = comsou["iers_name"][i]

    [count, dra_x_sim, ddec_x_sim,
     dra_k_sim, ddec_k_sim,
     dra_ka_sim, ddec_ka_sim,
     dra_g_sim, ddec_g_sim] = np.genfromtxt("../data/sim-pos/%s.pos" % souname,
                                            unpack=True)

    # Simulated offset
    # K - X
    dra_k_x_sim = dra_k_sim - dra_x_sim
    ddec_k_x_sim = ddec_k_sim - ddec_x_sim

    # Ka - X
    dra_ka_x_sim = dra_ka_sim - dra_x_sim
    ddec_ka_x_sim = ddec_ka_sim - ddec_x_sim

    # Gaia - X
    dra_g_x_sim = dra_g_sim - dra_x_sim
    ddec_g_x_sim = ddec_g_sim - ddec_x_sim

    # Recompute the offset
    # K - X
    dra_k_x_new = dra_k_x_sim + comsou["dra_k_ccl2"][i]
    ddec_k_x_new = ddec_k_x_sim + comsou["ddec_k_ccl2"][i]
    ang_sep_k_x_new = np.sqrt(dra_k_x_new**2 + ddec_k_x_new**2)
    pa_k_x_new = pa_calc(dra_k_x_new, ddec_k_x_new)

    # Ka - X
    dra_ka_x_new = dra_ka_x_sim + comsou["dra_ka_ccl2"][i]
    ddec_ka_x_new = ddec_ka_x_sim + comsou["ddec_ka_ccl2"][i]
    ang_sep_ka_x_new = np.sqrt(dra_ka_x_new**2 + ddec_ka_x_new**2)
    pa_ka_x_new = pa_calc(dra_ka_x_new, ddec_ka_x_new)

    # Gaia - X
    dra_g_x_new = dra_g_x_sim + comsou["dra_g_ccl2"][i]
    ddec_g_x_new = ddec_g_x_sim + comsou["ddec_g_ccl2"][i]
    ang_sep_g_x_new = np.sqrt(dra_g_x_new**2 + ddec_g_x_new**2)
    pa_g_x_new = pa_calc(dra_g_x_new, ddec_g_x_new)

    # Estimate the offset and PA
    ang_sep_k_sim, ang_sep_kmin, ang_sep_kmax = search_for_rho_peak(ang_sep_k_x_new)
    ang_sep_ka_sim, ang_sep_kamin, ang_sep_kamax = search_for_rho_peak(ang_sep_ka_x_new)
    ang_sep_g_sim, ang_sep_gmin, ang_sep_gmax = search_for_rho_peak(ang_sep_g_x_new)

    pak = comsou["pa_k_ccl2"][i]
    paka = comsou["pa_ka_ccl2"][i]
    pag = comsou["pa_g_ccl2"][i]

    pa_k_sim, pa_kmin, pa_kmax = search_for_pa_peak(pa_k_x_new, pak)
    pa_ka_sim, pa_kamin, pa_kamax = search_for_pa_peak(pa_ka_x_new, paka)
    pa_g_sim, pa_gmin, pa_gmax = search_for_pa_peak(pa_g_x_new, pag)

    # Write the results
    with open("../data/sim-offset/{:s}.offset".format(souname), "w") as fdata:
        print("# Simulated position offsets at X, K, Ka, and Gaia", file=fdata)
        print("# Columns  Meaning", file=fdata)
        print("# Simulated position offsets at X, K, Ka, and Gaia", file=fdata)
        print("# Columns  Meaning", file=fdata)
        print("#    0     Count", file=fdata)
        print("#    1     K-X offset (mas)", file=fdata)
        print("#    2     PA of K-X offset (deg)", file=fdata)
        print("#    3     Ka-X offset (mas)", file=fdata)
        print("#    4     PA of Ka-X offset (deg)", file=fdata)
        print("#    5     Gaia-X offset (mas)", file=fdata)
        print("#    6     PA of Gaia-X offset (deg)", file=fdata)
        print("# Created at {}".format(time.asctime()), file=fdata)

        N = len(dra_ka_x_new)
        for j in range(N):
            fmt = "%8d  " + "%.3f  %.0f  " * 3
            print(fmt % ((j+1), ang_sep_k_x_new[j], pa_k_x_new[j],
                         ang_sep_ka_x_new[j], pa_ka_x_new[j],
                         ang_sep_g_x_new[j], pa_g_x_new[j]), file=fdata)

    print("complete {:d}/{:d}!".format(i+1, len(comsou)))

    return [ang_sep_k_sim, ang_sep_kmin, ang_sep_kmax,
            pa_k_sim, pa_kmin, pa_kmax,
            ang_sep_ka_sim, ang_sep_kamin, ang_sep_kamax,
            pa_ka_sim, pa_kamin, pa_kamax,
            ang_sep_g_sim, ang_sep_gmin, ang_sep_gmax,
            pa_g_sim, pa_gmin, pa_gmax]


# -----------------------------  MAINS ---------------------------------
def simulate_pos():
    """
    """

    comsou = Table.read("../data/common-source-position0.fits")

    for i in range(len(comsou)):
        # X-band
        ra_err, dec_err, ra_dec_corr = comsou["ra_err_x", "dec_err_x", "ra_dec_corr_x"][i]
        dra_x_sim, ddec_x_sim = simulate_from_covariance(ra_err, dec_err, ra_dec_corr)

        # K-band
        ra_err, dec_err, ra_dec_corr = comsou["ra_err_k", "dec_err_k", "ra_dec_corr_k"][0]
        dra_k_sim, ddec_k_sim = simulate_from_covariance(ra_err, dec_err, ra_dec_corr)

        # Ka-band
        ra_err, dec_err, ra_dec_corr = comsou["ra_err_ka", "dec_err_ka", "ra_dec_corr_ka"][0]
        dra_ka_sim, ddec_ka_sim = simulate_from_covariance(ra_err, dec_err, ra_dec_corr)

        # Gaia
        ra_err, dec_err, ra_dec_corr = comsou["ra_err_g", "dec_err_g", "ra_dec_corr_g"][0]
        dra_g_sim, ddec_g_sim = simulate_from_covariance(ra_err, dec_err, ra_dec_corr)

        # Source name
        sou = comsou["iers_name"][i]

        # Write the results
        with open("../data/sim-pos/{:s}.pos".format(sou), "w") as fdata:
            print("# Simulated positions at X, K, Ka, and Gaia", file=fdata)
            print("# Columns  Meaning", file=fdata)
            print("#    0     Count", file=fdata)
            print("#    1     RA offset at X-band (mas)", file=fdata)
            print("#    2     Decl. offset at X-band (mas)", file=fdata)
            print("#    3     RA offset at K-band (mas)", file=fdata)
            print("#    4     Decl. offset at K-band (mas)", file=fdata)
            print("#    5     RA offset at Ka-band (mas)", file=fdata)
            print("#    6     Decl. offset at Ka-band (mas)", file=fdata)
            print("#    7     RA offset at Gaia (mas)", file=fdata)
            print("#    8     Decl. offset at Gaia (mas)", file=fdata)
            print("# Created at {}".format(time.asctime()), file=fdata)

            for j in range(len(dra_x_sim)):
                print("{:8d}  {:+7.3f}  {:+7.3f}  {:+7.3f}  {:+7.3f}  "
                      "{:+7.3f}  {:+7.3f}  {:+7.3f}  {:+7.3f}".format(
                          (j+1),
                          dra_x_sim[j], ddec_x_sim[j], dra_k_sim[j], ddec_k_sim[j],
                          dra_ka_sim[j], ddec_ka_sim[j], dra_g_sim[j], ddec_g_sim[j]),
                      file=fdata)

            print("Simulation of {:s}: done!".format(sou))


# Recompute the offset
def recalc_offset():
    """
    """

    comsou = Table.read("../data/common-source-position.fits")
    N = len(comsou)

    with open("../logs/sim_pa_and_offset.log", "w") as fdata:
        # Add header information
        print("# Simulated position offsets at X, K, Ka, and Gaia", file=fdata)
        print("# Columns  Meaning", file=fdata)
        print("#    0     IERS name", file=fdata)
        print("#    1     K-X offset (mas)", file=fdata)
        print("#    2     down limit (mas)", file=fdata)
        print("#    3     up limit (mas)", file=fdata)
        print("#    4     PA of K-X offset (deg)", file=fdata)
        print("#    5     down limit (deg)", file=fdata)
        print("#    6     up limit (deg)", file=fdata)
        print("#    7     Ka-X offset (mas)", file=fdata)
        print("#    8     down limit (mas)", file=fdata)
        print("#    9     up limit (mas)", file=fdata)
        print("#   10     PA of Ka-X offset (deg)", file=fdata)
        print("#   11     down limit (deg)", file=fdata)
        print("#   12     up limit (deg)", file=fdata)
        print("#   13     Gaia-X offset (mas)", file=fdata)
        print("#   14     down limit (mas)", file=fdata)
        print("#   15     up limit (mas)", file=fdata)
        print("#   16     PA of Gaia-X offset (deg)", file=fdata)
        print("#   17     down limit (deg)", file=fdata)
        print("#   18     up limit (deg)", file=fdata)
        print("# Created at {}".format(time.asctime()), file=fdata)

        for i in range(N):
            souname = comsou["iers_name"][i]
            data = calc_sim_offset(comsou, i)
            fmt = "%s  "+"%.3f  %.3f  %.3f  %.0f  %.0f  %.0f  " * 3
            print(fmt % (souname, *data), file=fdata)


def sim_hist():
    """
    """

    comsou = Table.read("../data/com-sou-list.txt", format="ascii")
    N = len(comsou)

    # Concate all the data
    souname = comsou["iers_name"][0]
    tab = Table.read("../data/sim-offset/{:s}.offset".format(souname), format="ascii",
                     names=["Count", "rho_k", "pa_k", "rho_ka", "pa_ka",
                            "rho_g", "pa_g"])
    print("Read data of ", souname)

    for num in range(1, N):
        souname = comsou["iers_name"][num]
        print("Read data of ", souname)
        tmp = Table.read("../data/sim-offset/{:s}.offset".format(souname), format="ascii",
                         names=["Count", "rho_k", "pa_k", "rho_ka", "pa_ka",
                                "rho_g", "pa_g"])
        tab = vstack([tab, tmp])

    # iter_time = 2
    iter_time = 10000

    # Output files
    fpa_k = open("../data/pa_k_sim_hist.dat", "w")
    fpa_ka = open("../data/pa_ka_sim_hist.dat", "w")
    fpa_g = open("../data/pa_gaia_sim_hist.dat", "w")
    fpa_k_g = open("../data/pa_k_gaia_sim_hist.dat", "w")
    fpa_ka_g = open("../data/pa_ka_gaia_sim_hist.dat", "w")
    fpa_k_ka = open("../data/pa_k_ka_sim_hist.dat", "w")

    for iter in range(iter_time):
        mask = (tab["Count"] == iter + 1)
        tab1 = tab[mask]

        pa_ki = np.array(tab1["pa_k"])
        pa_kai = np.array(tab1["pa_ka"])
        pa_gi = np.array(tab1["pa_g"])

        pa_k_gi = wrap_angle(pa_ki - pa_gi)
        pa_ka_gi = wrap_angle(pa_kai - pa_gi)
        pa_k_kai = wrap_angle(pa_ki - pa_kai)

        numbin = 15
        bins_set = np.linspace(0, 360, numbin)

        # PA (K-X)
        count, binedge = np.histogram(pa_ki, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_k)

        # PA (Ka-X)
        count, binedge = np.histogram(pa_kai, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_ka)

        # PA (Gaia-X)
        count, binedge = np.histogram(pa_gi, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_g)

        bins_set = np.linspace(-180, 180, numbin)
        # PA (K-X) vs PA (Gaia-X)
        count, binedge = np.histogram(pa_k_gi, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_k_g)

        # PA (Ka-X) vs PA (Gaia-X)
        count, binedge = np.histogram(pa_ka_gi, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_ka_g)

        # PA (K-X) vs PA (Ka-X)
        count, binedge = np.histogram(pa_k_kai, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_k_ka)

        print("{:10d} time ends".format(iter+1))

    fpa_k.close()
    fpa_ka.close()
    fpa_g.close()
    fpa_k_g.close()
    fpa_k_ka.close()
    fpa_ka_g.close()


# simulate_pos()
recalc_offset()
sim_hist()
# --------------------------------- END --------------------------------
