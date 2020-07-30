#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: get_astrogeo_image.py
"""
Created on Sun Jul 26 14:47:40 2020

@author: Neo(liuniu@smail.nju.edu.cn)

Use the web crawler technique to get the image data for sources we are
interested in.
"""

import re
from astropy.table import Table
from bs4 import BeautifulSoup
import requests


# -----------------------------  FUNCTIONS -----------------------------
def get_sou_image_link(sou_name):
    """Get the download links of images in Astrogeo database given source name
    """

    # Astrogeo image database link for query the source image
    url = "http://astrogeo.org/cgi-bin/imdb_get_source.csh"

    data = {"source_name": sou_name}

    content = requests.get(url=url, params=data).content

    soup = BeautifulSoup(content, features="lxml-xml")

    url_list = []

    # Check if the image is available
    if len(soup):

        # For all images
        tag_list = soup.find_all("A", HREF=re.compile("map.fits"))

        # For X_band image only
        # soup.find_all("A", string="X_map_fits")

        for tag in tag_list:
            url_list.append(tag["HREF"])

    return url_list


def parse_url_link(url_link):
    """
    """

    # A sample for link is
    # /images/J0130+0842/J0130+0842_S_1995_07_15_yyk_map.fits

    data = url_link.split("/")[-1].split("_")

    # J2000 name
    j2000_name = data[0]

    # band
    band = data[1]

    # epoch
    epoch = "{}-{}-{}".format(*data[2:5])

    # Author
    author = data[5]

    return j2000_name, band, epoch, author


# Source list
sou_list = Table.read("../data/miss_sou_in_bvid.txt", format="ascii")

# Webpage of Astrogeo
main_web = "http://astrogeo.org"


# Output file
with open("../data/astrogeo_image_info.txt", "w") as fop:

    # Add header
    print("iers_name,J2000_name,epoch,band,link,author",
          file=fop)

    for sou_name in sou_list["iers_name"]:
        url_list = get_sou_image_link(sou_name)

        if len(url_list):
            for url_link in url_list:
                j2000_name, band, epoch, author = parse_url_link(url_link)
                print("{},{},{},{},{},{}".format(sou_name, j2000_name, epoch,
                                                 band, main_web+url_link,
                                                 author), file=fop)
        else:
            print("There is no image for", sou_name)


# --------------------------------- END --------------------------------
