#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: classify_source.py
"""
Created on Mon Jul 27 09:03:11 2020

@author: Neo(liuniu@smail.nju.edu.cn)

According to availability of the images data in the BVID and Astrogeo database
to classify the sources.

"""

from astropy.table import Table, join
import numpy as np


# -----------------------------  MAIN -----------------------------
# Table of 488 sources
com_sou_list = Table.read("../data/com-sou-list.txt", format="ascii")

# BVID SI data
si_x_tab = Table.read("../data/si_x_0527.fits")

# Sources with BIVD image
mask = si_x_tab["si_c2015"] != 0
si_x_A = si_x_tab[mask]

# Sources without BIVD image
mask = si_x_tab["si_c2015"] == 0
si_x_B = si_x_tab[mask]

# Astrogeo data
img_db_info = Table.read("../logs/astrogeo-image.info", format="ascii")

# Sources with Astrogeo image between 2014 and 2016
img_x_1 = img_db_info[img_db_info["2015_flag"] == 1]


# Sources with Astrogeo image but not among 2014-2016
img_x_2 = img_db_info[img_db_info["2015_flag"] == 2]

with open("../logs/sou-class.txt", "w") as fop:

    # Add header information
    print("iers_name,class", file=fop)

    for sou in com_sou_list["iers_name"]:

        if sou in si_x_A["iers_name"]:
            flag = "A"

        elif sou in si_x_B["iers_name"]:
            if sou in img_x_1["iers_name"]:
                flag = "B1"
            elif sou in img_x_2["iers_name"]:
                flag = "B2"
            else:
                flag = "B3"

        else:
            if sou in img_x_1["iers_name"]:
                flag = "C1"
            elif sou in img_x_2["iers_name"]:
                flag = "C2"
            else:
                flag = "C3"

        print("{},{}".format(sou, flag), file=fop)


# --------------------------------- END --------------------------------
