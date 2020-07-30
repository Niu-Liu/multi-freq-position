#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: fig_copy.py
"""
Created on Mon May 13 00:37:57 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import os

# -----------------------------  FUNCTIONS -----------------------------
figfiles = ["0552+398.eps",
            "0859-140.eps",
            "1157-215.eps",
            "1243-072.eps",
            "1655+077.eps",
            "1742-078.eps",
            "2134+004.eps",
            "0430+052.eps",
            "0723-008.eps",
            "0743-006.eps",
            "0827+243.eps",
            "0850+581.eps",
            "2251+158.eps",
            "0003+380.eps",
            "0112-017.eps",
            "0122-003.eps",
            "0133+476.eps",
            "0917+449.eps",
            "0953+254.eps",
            "1038+064.eps",
            "1751+288.eps",
            "2234+282.eps",
            "k-sx-scatter.eps",
            "xka-sx-scatter.eps",
            "gaia-sx-scatter.eps",
            "rho-hist.eps",
            "X-hist.eps",
            "pa-hist.eps",
            # "rho-com-vs-si.eps",
            # "X-com-vs-si.eps",
            "rho-g-mag.eps",
            "X-g-mag.eps",
            "rho-bp-rp.eps",
            "rho-z.eps",
            "rho-I1R.eps",
            "rho-I2R.eps",
            "rho-I3R.eps",
            "rho-si.eps",
            "pa-diff.eps",
            "jet-pa-com.eps",
            #             "sim-pa-offset.png",
            # "image/1157-215.u.2013_05_05.icn_color.png",
            #             "sim-rho-vs.eps",
            #             "sim-pa-vs.eps",
            #             "0838+456-sim-x.png",
            #             "0838+456-sim-k.png",
            #             "0838+456-sim-ka.png",
            #             "0838+456-sim-gaia.png",
            #             "0838+456-sim-k-x.png",
            #             "0838+456-sim-ka-x.png",
            #             "0838+456-sim-g-x.png",
            #             "0838+456-sim-k-x-oft.png",
            #             "0838+456-sim-ka-x-oft.png",
            #             "0838+456-sim-g-x-oft.png"
            ]

dir1 = "../plots/"
# dir2 = "../notes/20191110/figs/"
# dir2 = "../notes/20200211/figs/"
# dir2 = "../notes/20200224/figs/"
dir2 = "../notes/20200407/figs/"

print("Input directory: ", dir1)
print("Output directory: ", dir2)

# Remove old figures
os.system("rm {}*".format(dir2))

for fig in figfiles:
    os.system("cp {}{} {}".format(dir1, fig, dir2))
    print("Copy {}: done!".format(fig))


# --------------------------------- END --------------------------------
