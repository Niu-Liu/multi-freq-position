#!/usr/local/bin/bash
#########################################################################
# File Name: find_more_sou_in_astrogeo.sh
# Author: Neo
# mail: niu.liu@nju.edu.cn
# Created Time: Sun Feb 21 16:48:58 2021
#########################################################################

ipython get_astrogeo_image.py
ipython classify_source.py
ipython astrogeo_image_stat.py
