#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: corr-test.py
"""
Created on Tue Jul 28 16:25:44 2020

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table
import numpy as np
from scipy.stats import spearmanr, pearsonr, kendalltau


# -----------------------------  MAIN -----------------------------
def CorrTest(x, y):
    """Correlation test between x and y
    """

    # Calculate the correlation between x and y
    # Pearson
    pcor = pearsonr(x, y)
    r, rp = pcor

    # Spearman
    scor = spearmanr(x, y)
    rs, rsp = scor

    # Kendall
    kcor = kendalltau(x, y)
    tau, taup = kcor

    # Combine the result
    # cor = {}
    # cor["r"] = r
    # cor["rp"] = rp
    # cor["rs"] = rs
    # cor["rsp"] = rsp
    # cor["tau"] = tau
    # cor["taup"] = taup
    cor = [r, rp, rs, rsp, tau, taup]

    return cor


def R2OCorrTest(tab, keys, bin_array):
    """
    """

    # Construct an empty Table to store the data
    zero_array = np.zeros_like(bin_array) * 1.0

    # For rho vs G
    rho_corr = Table([zero_array, zero_array,
                      zero_array, zero_array, zero_array, zero_array,
                      zero_array, zero_array, zero_array, zero_array,
                      zero_array, zero_array, zero_array, zero_array,
                      zero_array, zero_array, zero_array, zero_array],
                     names=["r_sx", "rp_sx", "rs_sx", "rsp_sx", "tau_sx", "taup_sx",
                            "r_k", "rp_k", "rs_k", "rsp_k", "tau_k", "taup_k",
                            "r_ka", "rp_ka", "rs_ka", "rsp_ka", "tau_ka", "taup_ka"])

    # For X vs G
    x_corr = Table(rho_corr)

    # A loop to check the choice of bin size
    number = np.arange(len(tab))
    for i, binsize in enumerate(bin_array):
        nb_bin = np.trunc(number / binsize)

        # group the table
        tab_g = tab.group_by(nb_bin)
        tab_bin = tab_g.groups.aggregate(np.median)

        # 1. rho vs G
        corr_sx_rho = CorrTest(tab_bin[keys], tab_bin["ang_sep_sx_1"])
        corr_k_rho = CorrTest(tab_bin[keys], tab_bin["ang_sep_k_1"])
        corr_ka_rho = CorrTest(tab_bin[keys], tab_bin["ang_sep_ka_1"])

        # 1.1. SX-band
        # 1.1.1 Pearson test
        rho_corr["r_sx"][i], rho_corr["rp_sx"][i] = corr_sx_rho[:2]
        # 1.1.2 Spearman
        rho_corr["rs_sx"][i], rho_corr["rsp_sx"][i] = corr_sx_rho[2:4]
        # 1.1.3 Kendall
        rho_corr["tau_sx"][i], rho_corr["taup_sx"][i] = corr_sx_rho[4:]

        # 1.2. K-band
        # 1.2.1 Pearson test
        rho_corr["r_k"][i], rho_corr["rp_k"][i] = corr_k_rho[:2]
        # 1.2.2 Spearman
        rho_corr["rs_k"][i], rho_corr["rsp_k"][i] = corr_k_rho[2:4]
        # 1.2.3 Kendall
        rho_corr["tau_k"][i], rho_corr["taup_k"][i] = corr_k_rho[4:]

        # 1.3. K-band
        # 1.3.1 Pearson test
        rho_corr["r_ka"][i], rho_corr["rp_ka"][i] = corr_ka_rho[:2]
        # 1.3.2 Spearman
        rho_corr["rs_ka"][i], rho_corr["rsp_ka"][i] = corr_ka_rho[2:4]
        # 1.3.3 Kendall
        rho_corr["tau_ka"][i], rho_corr["taup_ka"][i] = corr_ka_rho[4:]

        # 2. X vs G
        corr_sx_x = CorrTest(tab_bin[keys], tab_bin["nor_sep_sx_1"])
        corr_k_x = CorrTest(tab_bin[keys], tab_bin["nor_sep_k_1"])
        corr_ka_x = CorrTest(tab_bin[keys], tab_bin["nor_sep_ka_1"])

        # 2.1. SX-band
        # 2.1.1 Pearson test
        x_corr["r_sx"][i], x_corr["rp_sx"][i] = corr_sx_x[:2]
        # 2.1.2 Spearman
        x_corr["rs_sx"][i], x_corr["rsp_sx"][i] = corr_sx_x[2:4]
        # 2.1.3 Kendall
        x_corr["tau_sx"][i], x_corr["taup_sx"][i] = corr_sx_x[4:]

        # 2.2. K-band
        # 2.2.1 Pearson test
        x_corr["r_k"][i], x_corr["rp_k"][i] = corr_k_x[:2]
        # 2.2.2 Spearman
        x_corr["rs_k"][i], x_corr["rsp_k"][i] = corr_k_x[2:4]
        # 2.2.3 Kendall
        x_corr["tau_k"][i], x_corr["taup_k"][i] = corr_k_x[4:]

        # 2.3. XKa-band
        # 2.3.1 Pearson test
        x_corr["r_ka"][i], x_corr["rp_ka"][i] = corr_ka_x[:2]
        # 2.3.2 Spearman
        x_corr["rs_ka"][i], x_corr["rsp_ka"][i] = corr_ka_x[2:4]
        # 2.3.3 Kendall
        x_corr["tau_ka"][i], x_corr["taup_ka"][i] = corr_ka_x[4:]

    return rho_corr, x_corr


if __name__ == "__main__":
    main()
# --------------------------------- END --------------------------------
