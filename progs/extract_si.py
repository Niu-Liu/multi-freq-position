#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: extract_si.py
"""
Created on Tue Jul 21 09:59:23 2020

@author: Neo(liuniu@smail.nju.edu.cn)

Extract the Structure index data from those downloaded from BVID

I use two statistics as the representation on the SI data:
1) the one made at the median epoch of all time
2) the one closest to J2015.0
"""


# -----------------------------  MAIN -----------------------------

from astropy.table import Table, join
import numpy as np


def si_data_handle(si_tab, sou_list, band="s"):
    """Handle the SI data
    """

    # Cross-match
    si_4_com_sou = join(si_tab, sou_list, keys="iers_name")

    # To see how many sources are given a SI measurement
    si_4_com_sou_g = si_4_com_sou.group_by("iers_name")

    com_sou_with_si = si_4_com_sou_g.groups.keys
    print(len(com_sou_with_si),
          "sources are given the SI seasurements at {}-band.".format(band))

    # Empty array to store data
    num_sou_with_si = len(com_sou_with_si)
    si_med = np.zeros(num_sou_with_si)
    epo_med = np.zeros(num_sou_with_si)
    si_c2015 = np.zeros(num_sou_with_si)
    epo_c2015 = np.zeros(num_sou_with_si)

    for i, group in enumerate(si_4_com_sou_g.groups):

        if len(group) == 1:
            # Median value
            si_med[i] = group["si"][0]
            epo_med[i] = group["year"][0]

            # value near J2015
            if np.fabs(group["year"][0] - 2015) <= 1:
                si_c2015[i] = group["si"][0]
                epo_c2015[i] = group["year"][0]
            else:
                si_c2015[i] = 0
                epo_c2015[i] = 0

        else:
            # Median value
            sort_ind = np.argsort(group["si"])
            si_sort = group["si"][sort_ind]
            epo_sort = group["year"][sort_ind]
            med_ind = int(len(group) / 2)
            si_med[i] = si_sort[med_ind]
            epo_med[i] = epo_sort[med_ind]

            # Value nearest J2015
            dif_2015 = np.fabs(group["year"] - 2015)
            # We need the SI made within 1-yr to J2015
            if np.min(dif_2015) <= 1:
                # ind = dif_2015.argmin()
                # si_c2015[i] = group["si"][ind]
                # epo_c2015[i] = group["year"][ind]

                # Or use the median value
                mask = (group["year"] < 2016) & (group["year"] > 2014)
                si_c2015[i] = np.median(group[mask]["si"])
                epo_c2015[i] = np.median(group[mask]["year"])

            else:
                si_c2015[i] = 0
                epo_c2015[i] = 0

    si_tab = Table([com_sou_with_si["iers_name"],
                    si_med, epo_med, si_c2015, epo_c2015],
                   names=["iers_name", "si_med", "epo_med", "si_c2015", "epo_c2015"])

    return si_tab


# SI data
# S-band image
si_s_table = Table.read("/Users/Neo/Astronomy/data/bvid/SI_S_20200527.fits")
si_s_table.rename_columns(["Source", "Date", "Flux (Jy)", "Structure index", "Compactness"],
                          ["iers_name", "year", "flux", "si", "compactness"])

# X-band image
si_x_table = Table.read("/Users/Neo/Astronomy/data/bvid/SI_X_20200527.fits")
si_x_table.rename_columns(["Source", "Date", "Flux (Jy)", "Structure index", "Compactness"],
                          ["iers_name", "year", "flux", "si", "compactness"])

# K-band image
si_k_table = Table.read("/Users/Neo/Astronomy/data/bvid/SI_K_20200527.fits")
si_k_table.rename_columns(["Source", "Date", "Flux (Jy)", "Structure index", "Compactness"],
                          ["iers_name", "year", "flux", "si", "compactness"])

# Q-band image
si_q_table = Table.read("/Users/Neo/Astronomy/data/bvid/SI_Q_20200527.fits")
si_q_table.rename_columns(["Source", "Date", "Flux (Jy)", "Structure index", "Compactness"],
                          ["iers_name", "year", "flux", "si", "compactness"])

# Table of 488 sources
com_sou_list = Table.read("../data/com-sou-list.txt", format="ascii")

# Work on SI data
si_s_tab = si_data_handle(si_s_table, com_sou_list, band="S")
si_x_tab = si_data_handle(si_x_table, com_sou_list, band="X")
si_k_tab = si_data_handle(si_k_table, com_sou_list, band="K")
si_q_tab = si_data_handle(si_q_table, com_sou_list, band="Q")

# Save those data.
si_s_tab.write("../data/si_s_0527.fits", overwrite=True)
si_x_tab.write("../data/si_x_0527.fits", overwrite=True)
si_k_tab.write("../data/si_k_0527.fits", overwrite=True)
si_q_tab.write("../data/si_q_0527.fits", overwrite=True)

# --------------------------------- END --------------------------------
