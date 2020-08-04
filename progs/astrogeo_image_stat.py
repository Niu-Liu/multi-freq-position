#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: astrogeo-image-stat.py
"""
Created on Mon Jul 27 09:26:13 2020

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table, setdiff
from astropy.time import Time
import numpy as np

# -----------------------------  MAIN -----------------------------
# Table of 488 sources
# com_sou_list = Table.read("../data/com-sou-list.txt", format="ascii")

img_db_info = Table.read("../data/astrogeo_image_info.txt", format="ascii")

# Keep images at X-band
mask = img_db_info["band"] == "X"
x_img_db_info = img_db_info[mask]

# Epoch from string to astropy.time.Time object
epoch = Time(x_img_db_info["epoch"], scale="utc")
x_img_db_info["epoch"] = epoch.jyear

# Group by source name
x_img_db_info_g = x_img_db_info.group_by("iers_name")
# sou_list = x_img_db_info_g.groups.keys

# Output file to store data
with open("../logs/astrogeo-image.info", "w") as fop:

    # Add header information
    print("iers_name,2015_flag", file=fop)

    for i, group in enumerate(x_img_db_info_g.groups):

        # if len(group["epoch"]) and group["epoch"][0]:
        epo_dif = np.fabs(group["epoch"] - 2015)
        epo_dif_min = np.min(epo_dif)

        if epo_dif_min <= 1:
            flag = 1
        else:
            flag = 2

        # else:
            # flag = 3

        print("{},{:.0f}".format(group["iers_name"][0], flag), file=fop)

# print X-band image
with open("../logs/astrogeo-x-image.dat", "w") as fop:
    for link in x_img_db_info["link"]:
        print(link, file=fop)


# --------------------------------- END --------------------------------
