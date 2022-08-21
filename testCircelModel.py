#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 08:12:36 2021

@author: svlesovoi
"""

import numpy as NP;
import pylab as PL;
from astropy.io import fits
from BadaryRAO import BadaryRAO
from sunpy import coordinates
from matplotlib.ticker import (MultipleLocator)
import skimage
from skimage import filters
from skimage import transform
from skimage import measure
from skimage.feature import canny
from skimage.transform import warp, AffineTransform
from matplotlib.patches import Ellipse, Circle

L = 1024
radius = L/8
x0 = L/2-L/5
y0 = L/2-L/3
l0 = L/2
image = NP.zeros((L,L))

for i in range(L):
    for j in range(L):
        if (NP.sqrt((i - x0)**2 + (j - y0)**2) < radius):
            image[i,j] = 1

sunCir = skimage.measure.CircleModel()
contours = measure.find_contours(image, 0.5)
contourLength = []
for n, contour in enumerate(contours):
    contourLength.append(len(contour))
contourLength = NP.array(contourLength)

maxInd = NP.argmax(contourLength)
sunCir.estimate(contours[maxInd])
i_cx0, i_cy0, i_cR = sunCir.params
iImage = warp(image,AffineTransform(translation=(l0 - i_cy0,l0 - i_cx0)).inverse)

PL.figure()
PL.imshow(image)

PL.figure()
PL.imshow(iImage)
    

