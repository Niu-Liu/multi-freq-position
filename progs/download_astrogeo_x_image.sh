#!/usr/local/bin/bash
#########################################################################
# File Name: download_astrogeo_x_image.sh
# Author: Neo
# mail: liuniu@smail.nju.edu.cn
# Created Time: Tue Aug  4 16:27:48 2020
#########################################################################

dest_dir="../data/astrogeo-x-image"
mkdir ${dest_dir} 
cd ${dest_dir}

for link in `cat ../../logs/astrogeo-x-image.dat`
do
    wget ${link}
done
