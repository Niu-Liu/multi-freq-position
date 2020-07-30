#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: calculate-radio-to-optical-distance-with-ICRF3-SX.py
"""
Created on Mon Jul 27 19:03:23 2020

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table, join, Column
import astropy.units as u
import numpy as np
# My modules
# from my_progs.catalog.pos_diff import nor_sep_calc
from my_progs.catalog import read_icrf, read_gaia, pos_diff
from calc_pa import pa_calc

nor_sep_calc = pos_diff.nor_sep_calc

# -----------------------------  MAIN -----------------------------
# Offsets wrt. SX position
k2x = Table.read("../data/icrf3_k_sx_offset.fits")
ka2x = Table.read("../data/icrf3_xka_sx_offset.fits")
g2x = Table.read("../data/gaia_icrf3_sx_offset.fits")

# Keep columns needed in the following calculation
# K - SX
k2x.keep_columns(["iers_name", "ra_err_k", "dec_err_k",
                  "ra", "dec", "ra_err_x", "dec_err_x",
                  "dra_k", "ddec_k", "dra_k_ccl2", "ddec_k_ccl2"])
# Resort the table
new_order = ["iers_name", "ra", "dec",
             "ra_err_x", "dec_err_x", "ra_err_k", "dec_err_k",
             "dra_k", "ddec_k", "dra_k_ccl2", "ddec_k_ccl2"]
k2x = k2x[new_order]
k2x.rename_columns(["dra_k_ccl2", "ddec_k_ccl2"], ["dra_k_1", "ddec_k_1"])

# XKa - SX
ka2x.keep_columns(["iers_name", "ra_err_ka", "dec_err_ka",
                   "dra_ka", "ddec_ka", "dra_ka_ccl2", "ddec_ka_ccl2"])
# Resort the table
new_order = ["iers_name", "ra_err_ka", "dec_err_ka",
             "dra_ka", "ddec_ka", "dra_ka_ccl2", "ddec_ka_ccl2"]
ka2x = ka2x[new_order]
ka2x.rename_columns(["dra_ka_ccl2", "ddec_ka_ccl2"], ["dra_ka_1", "ddec_ka_1"])


# Gaia - SX
g2x.keep_columns(["iers_name", "ra_err_g", "dec_err_g",
                  "dra_g", "ddec_g", "ang_sep_g", "pa_g", "nor_sep_g",
                  "dra_g_ccl2", "ddec_g_ccl2",
                  "ang_sep_g_ccl2", "pa_g_ccl2", "nor_sep_g_ccl2"])
# Resort the table
new_order = ["iers_name", "ra_err_g", "dec_err_g",
             "dra_g", "ddec_g", "ang_sep_g", "pa_g", "nor_sep_g",
             "dra_g_ccl2", "ddec_g_ccl2",
             "ang_sep_g_ccl2", "pa_g_ccl2", "nor_sep_g_ccl2"]
g2x = g2x[new_order]
g2x.rename_columns(["dra_g_ccl2", "ddec_g_ccl2",
                    "ang_sep_g_ccl2", "pa_g_ccl2", "nor_sep_g_ccl2"],
                   ["dra_g_1", "ddec_g_1", "ang_sep_g_1", "pa_g_1", "nor_sep_g_1"])

# We still need the correlation between RA and Decl.
# Read catalogs
icrf3sx = read_icrf.read_icrf3(wv="sx")
icrf3k = read_icrf.read_icrf3(wv="k")
icrf3xka = read_icrf.read_icrf3(wv="xka")
gaiadr2 = read_gaia.read_dr2_iers()

icrf3sx.keep_columns(["iers_name", "ra_dec_corr"])
icrf3sx.rename_column("ra_dec_corr", "ra_dec_corr_x")
icrf3k.keep_columns(["iers_name", "ra_dec_corr"])
icrf3k.rename_column("ra_dec_corr", "ra_dec_corr_k")
icrf3xka.keep_columns(["iers_name", "ra_dec_corr"])
icrf3xka.rename_column("ra_dec_corr", "ra_dec_corr_ka")
gaiadr2.keep_columns(["iers_name", "ra_dec_corr"])
gaiadr2.rename_column("ra_dec_corr", "ra_dec_corr_g")

# Cross-match amongs tables
r2x = join(k2x, ka2x, keys="iers_name")
r2x = join(r2x, g2x, keys="iers_name")
tab1 = join(icrf3sx, icrf3k, keys="iers_name")
tab2 = join(icrf3xka, gaiadr2, keys="iers_name")
tab3 = join(tab1, tab2, keys="iers_name")
r2x = join(r2x, tab3, keys="iers_name")

# Calculate the radio-to-optical offset
# 1. K - Gaia
# 1.1 pre-fit
dra_k_g = r2x["dra_g"] - r2x["dra_k"]
ddec_k_g = r2x["ddec_g"] - r2x["ddec_k"]
dra_err_k_g = np.sqrt(r2x["ra_err_k"]**2 + r2x["ra_err_g"]**2)
ddec_err_k_g = np.sqrt(r2x["dec_err_k"]**2 + r2x["dec_err_g"]**2)
cov_k_g = r2x["ra_err_k"] * r2x["dec_err_k"] * r2x["ra_dec_corr_k"] + \
    r2x["ra_err_g"] * r2x["dec_err_g"] * r2x["ra_dec_corr_g"]
cor_k_g = cov_k_g / dra_err_k_g / ddec_err_k_g

ang_sep_k_g, Xa_k_g, Xd_k_g, X_k_g = nor_sep_calc(
    dra_k_g, dra_err_k_g, ddec_k_g, ddec_err_k_g, cor_k_g)
PA_k_g = pa_calc(dra_k_g, ddec_k_g)
PA_k_g = Column(PA_k_g, unit=u.deg)

# 1.2 post-fit
dra_k_g_1 = r2x["dra_g_1"] - r2x["dra_k_1"]
ddec_k_g_1 = r2x["ddec_g_1"] - r2x["ddec_k_1"]
# The covariance between RA and decl. of post-fit is assumed to be same with pre-fit.

ang_sep_k_g_1, Xa_k_g_1, Xd_k_g_1, X_k_g_1 = nor_sep_calc(
    dra_k_g_1, dra_err_k_g, ddec_k_g_1, ddec_err_k_g, cor_k_g)
PA_k_g_1 = pa_calc(dra_k_g_1, ddec_k_g_1)
PA_k_g_1 = Column(PA_k_g_1, unit=u.deg)

# 2. XKa - Gaia
# 2.1 pre-fit
dra_ka_g = r2x["dra_g"] - r2x["dra_ka"]
ddec_ka_g = r2x["ddec_g"] - r2x["ddec_ka"]
dra_err_ka_g = np.sqrt(r2x["ra_err_ka"]**2 + r2x["ra_err_g"]**2)
ddec_err_ka_g = np.sqrt(r2x["dec_err_ka"]**2 + r2x["dec_err_g"]**2)
cov_ka_g = r2x["ra_err_ka"] * r2x["dec_err_ka"] * r2x["ra_dec_corr_ka"] + \
    r2x["ra_err_g"] * r2x["dec_err_g"] * r2x["ra_dec_corr_g"]
cor_ka_g = cov_ka_g / dra_err_ka_g / ddec_err_ka_g

ang_sep_ka_g, Xa_ka_g, Xd_ka_g, X_ka_g = nor_sep_calc(
    dra_ka_g, dra_err_ka_g, ddec_ka_g, ddec_err_ka_g, cor_ka_g)
PA_ka_g = pa_calc(dra_ka_g, ddec_ka_g)
PA_ka_g = Column(PA_ka_g, unit=u.deg)

# 2.2 post-fit
dra_ka_g_1 = r2x["dra_g_1"] - r2x["dra_ka_1"]
ddec_ka_g_1 = r2x["ddec_g_1"] - r2x["ddec_ka_1"]
# The covariance between RA and decl. of post-fit is assumed to be same with pre-fit.
# This is the same with K-band

ang_sep_ka_g_1, Xa_ka_g_1, Xd_ka_g_1, X_ka_g_1 = nor_sep_calc(
    dra_ka_g_1, dra_err_ka_g, ddec_ka_g_1, ddec_err_ka_g, cor_ka_g)
PA_ka_g_1 = pa_calc(dra_ka_g_1, ddec_ka_g_1)
PA_ka_g_1 = Column(PA_ka_g_1, unit=u.deg)

# Construct a new table
r2o = r2x["iers_name", "ra", "dec",
          "dra_g", "ddec_g", "ang_sep_g", "pa_g", "nor_sep_g",
          "dra_g_1", "ddec_g_1", "ang_sep_g_1", "pa_g_1", "nor_sep_g_1"]
r2o.rename_columns(["dra_g", "ddec_g", "ang_sep_g", "pa_g", "nor_sep_g"],
                   ["dra_sx", "ddec_sx", "ang_sep_sx", "pa_sx", "nor_sep_sx"])
r2o.rename_columns(["dra_g_1", "ddec_g_1", "ang_sep_g_1", "pa_g_1", "nor_sep_g_1"],
                   ["dra_sx_1", "ddec_sx_1", "ang_sep_sx_1", "pa_sx_1", "nor_sep_sx_1"])
r2o.add_columns([dra_k_g, ddec_k_g, ang_sep_k_g, PA_k_g, X_k_g,
                 dra_ka_g, ddec_ka_g, ang_sep_ka_g, PA_ka_g, X_ka_g],
                names=["dra_k", "ddec_k", "ang_sep_k", "pa_k", "nor_sep_k",
                       "dra_ka", "ddec_ka", "ang_sep_ka", "pa_ka", "nor_sep_ka"])
r2o.add_columns([dra_k_g_1, ddec_k_g_1, ang_sep_k_g_1, PA_k_g_1, X_k_g_1,
                 dra_ka_g_1, ddec_ka_g_1, ang_sep_ka_g_1, PA_ka_g_1, X_ka_g_1],
                names=["dra_k_1", "ddec_k_1", "ang_sep_k_1", "pa_k_1", "nor_sep_k_1",
                       "dra_ka_1", "ddec_ka_1", "ang_sep_ka_1", "pa_ka_1", "nor_sep_ka_1"])

# Store the data
r2o.write("../data/multiwav-offset-in-SX-frame.fits", overwrite=True)
# --------------------------------- END --------------------------------
