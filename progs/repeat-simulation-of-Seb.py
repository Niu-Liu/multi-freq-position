#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: repeat-simulation-of-Seb.py
"""
Created on Tue Feb 25 10:22:02 2020

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table, join
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from numpy import sqrt, cos
from random import shuffle


from my_progs.catalog.read_gaia import read_dr2_iers
from my_progs.catalog.read_icrf import read_icrf3
from my_progs.catalog.vsh_deg2_cor import vsh_func_calc02
from calc_pa import pa_calc, wrap_angle

# -----------------------------  FUNCTIONS -----------------------------
method_1 = True
method_2 = False

gdr2 = read_dr2_iers()

icrf3x = read_icrf3(wv="sx")
icrf3k = read_icrf3(wv="k")
icrf3ka = read_icrf3(wv="xka")

tmp1 = join(icrf3x, icrf3k, keys="iers_name", table_names=["x", "k"])
tmp2 = join(icrf3ka, gdr2, keys="iers_name", table_names=["ka", "g"])

comsou = join(tmp1, tmp2, keys="iers_name")

comsou.keep_columns(["ra_x", "dec_x", "ra_k", "dec_k",
                     "ra_ka", "dec_ka", "ra_g", "dec_g"])

fpa_k = open("../data/sim-as-seb/pa_k_sim_hist_1.dat", "w")
fpa_ka = open("../data/sim-as-seb/pa_ka_sim_hist_1.dat", "w")
fpa_g = open("../data/sim-as-seb/pa_gaia_sim_hist_1.dat", "w")

fpa_k_g = open("../data/sim-as-seb/pa_k_gaia_sim_hist_1.dat", "w")
fpa_ka_g = open("../data/sim-as-seb/pa_ka_gaia_sim_hist_1.dat", "w")
fpa_k_ka = open("../data/sim-as-seb/pa_k_ka_sim_hist_1.dat", "w")

#iter_time = 1000
iter_time = 10000
numbin = 15

if method_1:
    for iteri in range(iter_time):
        # Histogram
        bins_set = np.linspace(0, 360, numbin)

        # Simulate K - X
        comsou1 = comsou["ra_x", "dec_x", "ra_k", "dec_k"]
        # Mess the X-band data
        shuffle(comsou1["ra_x"])
        shuffle(comsou1["dec_x"])

        # Mess the K-band data
        shuffle(comsou1["ra_k"])
        shuffle(comsou1["dec_k"])

        # Calculate PA
        coord_x = SkyCoord(comsou1["ra_x"], comsou1["dec_x"], frame="icrs")
        coord_k = SkyCoord(comsou1["ra_k"], comsou1["dec_k"], frame="icrs")
        pa = coord_x.position_angle(coord_k)
        pa_k_x = pa.to(u.deg).value

        # PA (K-X)
        count, binedge = np.histogram(pa_k_x, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_k)

        # Simulate Ka - X
        comsou2 = comsou["ra_x", "dec_x", "ra_ka", "dec_ka"]
        # Mess the X-band data
        shuffle(comsou2["ra_x"])
        shuffle(comsou2["dec_x"])

        # Mess the Ka-band data
        shuffle(comsou2["ra_ka"])
        shuffle(comsou2["dec_ka"])

        # Calculate PA
        coord_x = SkyCoord(comsou2["ra_x"], comsou2["dec_x"], frame="icrs")
        coord_ka = SkyCoord(comsou2["ra_ka"], comsou2["dec_ka"], frame="icrs")
        pa = coord_x.position_angle(coord_ka)
        pa_ka_x = pa.to(u.deg).value

        # PA (Ka-X)
        count, binedge = np.histogram(pa_ka_x, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_ka)

        # Simulate Gaia - X
        comsou3 = comsou["ra_x", "dec_x", "ra_g", "dec_g"]
        # Mess the X-band data
        shuffle(comsou3["ra_x"])
        shuffle(comsou3["dec_x"])

        # Mess the Gaia data
        shuffle(comsou3["ra_g"])
        shuffle(comsou3["dec_g"])

        # Calculate PA
        coord_x = SkyCoord(comsou3["ra_x"], comsou3["dec_x"], frame="icrs")
        coord_g = SkyCoord(comsou3["ra_g"], comsou3["dec_g"], frame="icrs")
        pa = coord_x.position_angle(coord_g)
        pa_g_x = pa.to(u.deg).value

        # PA (Gaia-X)
        count, binedge = np.histogram(pa_g_x, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_g)

        bins_set = np.linspace(-180, 180, numbin)

        # PA offset
        pa_k_g = wrap_angle(pa_k_x - pa_g_x)
        pa_k_ka = wrap_angle(pa_k_x - pa_ka_x)
        pa_ka_g = wrap_angle(pa_ka_x - pa_g_x)

        # PA (K-X) vs PA (Gaia-X)
        count, binedge = np.histogram(pa_k_g, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_k_g)

        # PA (Ka-X) vs PA (Gaia-X)
        count, binedge = np.histogram(pa_ka_g, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_ka_g)

        # PA (K-X) vs PA (Ka-X)
        count, binedge = np.histogram(pa_k_ka, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_k_ka)

        print("{:4d}th iteration: done!".format(iteri+1))

if method_2:

    # offset
    coord_x = SkyCoord(comsou["ra_x"], comsou["dec_x"], frame="icrs")
    coord_k = SkyCoord(comsou["ra_k"], comsou["dec_k"], frame="icrs")
    coord_ka = SkyCoord(comsou["ra_ka"], comsou["dec_ka"], frame="icrs")
    coord_g = SkyCoord(comsou["ra_g"], comsou["dec_g"], frame="icrs")

    # Offset
    dra_k, ddec_k = coord_x.spherical_offsets_to(coord_k)
    dra_ka, ddec_ka = coord_x.spherical_offsets_to(coord_ka)
    dra_g, ddec_g = coord_x.spherical_offsets_to(coord_g)

    cen = SkyCoord(0*u.deg, 0*u.deg, frame="icrs")

    for iteri in range(iter_time):
        # Mess position
        dra_k1 = shuffle(dra_k)
        ddec_k1 = shuffle(ddec_k)

        dra_ka1 = shuffle(dra_ka)
        ddec_ka1 = shuffle(ddec_ka)

        dra_g1 = shuffle(dra_g)
        ddec_g1 = shuffle(ddec_g)

        # Calculate PA
        oft_k = SkyCoord(dra_k, ddec_k, frame="icrs")
        oft_ka = SkyCoord(dra_ka, ddec_ka, frame="icrs")
        oft_g = SkyCoord(dra_g, ddec_g, frame="icrs")

        # K - X
        pa = cen.position_angle(oft_k)
        pa_k_x = pa.to(u.deg).value

        # Ka - X
        pa = cen.position_angle(oft_ka)
        pa_ka_x = pa.to(u.deg).value

        # Gaia - X
        pa = cen.position_angle(oft_g)
        pa_g_x = pa.to(u.deg).value

        # PA offset
        pa_k_g = wrap_angle(pa_k_x - pa_g_x)
        pa_k_ka = wrap_angle(pa_k_x - pa_ka_x)
        pa_ka_g = wrap_angle(pa_ka_x - pa_g_x)

        # Histogram
        bins_set = np.linspace(0, 360, numbin)

        # PA (K-X)
        count, binedge = np.histogram(pa_k_x, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_k)

        # PA (Ka-X)
        count, binedge = np.histogram(pa_ka_x, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_ka)

        # PA (Gaia-X)
        count, binedge = np.histogram(pa_g_x, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_g)

        bins_set = np.linspace(-180, 180, numbin)
        # PA (K-X) vs PA (Gaia-X)
        count, binedge = np.histogram(pa_k_g, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_k_g)

        # PA (Ka-X) vs PA (Gaia-X)
        count, binedge = np.histogram(pa_ka_g, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))
        print("".join(line), file=fpa_ka_g)

        # PA (K-X) vs PA (Ka-X)
        count, binedge = np.histogram(pa_k_ka, bins_set)
        line = []
        for counti in count:
            line.append("{:3d}  ".format(counti))

        print("".join(line), file=fpa_k_ka)

        print("{:4d}th iteration: done!".format(iteri+1))


fpa_k.close()
fpa_ka.close()
fpa_g.close()
fpa_k_g.close()
fpa_k_ka.close()
fpa_ka_g.close()


# --------------------------------- END --------------------------------
