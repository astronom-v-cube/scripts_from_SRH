#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 09:11:28 2020

@author: svlesovoi
"""

import os, fnmatch;
import numpy as NP;
import pylab as PL;
from astropy.io import fits
from BadaryRAO import BadaryRAO
from sunpy import coordinates
from matplotlib.ticker import (MultipleLocator)
import skimage
from skimage import filters
from skimage import transform
from matplotlib.patches import Ellipse

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename))
    result.sort()
    return result

srhFitsPath = '/home/svlesovoi/Documents/Python Scripts/srhFits'
fitNames = findFits(srhFitsPath, '*.fits')

theta = NP.deg2rad(NP.linspace(0,360,360))
sunEll = skimage.measure.EllipseModel()
radiusColors = PL.get_cmap('rainbow')
freqColors = radiusColors(NP.linspace(1.,0.,len(fitNames)))
fig, pl = PL.subplots()
pl.set_xlim(-20,20)
pl.set_ylim(-20,20)
freq = 0

eqProfiles = []
polProfiles = []
srhImages = []
srhEdges = []

for fitName in fitNames:
    srhFits = fits.open(fitName)
    srh_x_size = srhFits[0].header['NAXIS1']
    srh_y_size = srhFits[0].header['NAXIS2']
    srh_x_delt = srhFits[0].header['CDELT1'] / 60 * 2
    srh_y_delt = srhFits[0].header['CDELT2'] / 60 * 2
#    srhFitsData = NP.flipud(srhFits[0].data)
    srhFitsData = srhFits[0].data
    srhImages.append(srhFitsData)
    polProfiles.append(srhFitsData[:,int(srh_y_size//2)])
    eqProfiles.append(srhFitsData[int(srh_x_size//2),:])
    x0 = int(srh_x_size//2)
    y0 = int(srh_y_size//2)
    dx = int(srh_x_size//4)
    dy = int(srh_y_size//4)
    srhFitsDataMean = NP.mean(srhFitsData[x0 - dx:x0 + dx, y0 - dy:y0 + dy])
    srhFitsEdge = filters.sobel(NP.clip(srhFitsData,srhFitsDataMean/2 - srhFitsDataMean/100,srhFitsDataMean/2 + srhFitsDataMean/100))
    srhEdges.append(srhFitsEdge)
    ellInd = NP.where(srhFitsEdge > .5*srhFitsEdge.max())
    ellXY = NP.zeros((len(ellInd[1]),2))
    ellXY[:,0] = ellInd[1]
    ellXY[:,1] = ellInd[0]
    ellXY[:,0] = (ellXY[:,0] - x0) * srhFits[0].header['CDELT1'] / 60 * 2
    ellXY[:,1] = (ellXY[:,1] - y0) * srhFits[0].header['CDELT2'] / 60 * 2
    sunEll.estimate(ellXY)
    xc, yc, a, b, theta = sunEll.params
    transVector = transform.AffineTransform(translation=(-xc / srh_y_delt, -yc / srh_x_delt))
    srhImages[-1] = transform.warp(srhImages[-1], transVector.inverse)
    srhEdges[-1] = transform.warp(srhEdges[-1], transVector.inverse)
    print(xc, yc, a, b)
    ell_patch = Ellipse((0, 0), 2*a, 2*b, theta*180/NP.pi, edgecolor=freqColors[freq], facecolor='none', linewidth=1.)
    pl.add_patch(ell_patch)
    freq += 1
pl.grid()
#pl.imshow(srhEdges[0], extent=[-srh_x_size/2*srh_x_delt,srh_x_size/2*srh_x_delt,srh_y_size/2*srh_x_delt,-srh_y_size/2*srh_x_delt])
#pl.imshow(srhImages[3] + srhEdges[3]*100, extent=[-srh_x_size/2*srh_x_delt,srh_x_size/2*srh_x_delt,srh_y_size/2*srh_x_delt,-srh_y_size/2*srh_x_delt])
#pl.imshow(srhEdges[0]+srhEdges[1]+srhEdges[2]+srhEdges[3], extent=[-srh_x_size/2*srh_x_delt,srh_x_size/2*srh_x_delt,srh_y_size/2*srh_x_delt,-srh_y_size/2*srh_x_delt])
pl.imshow(srhImages[0]+srhImages[1]*2+srhImages[2]*3+srhImages[3]*4, extent=[-srh_x_size/2*srh_x_delt,srh_x_size/2*srh_x_delt,srh_y_size/2*srh_x_delt,-srh_y_size/2*srh_x_delt])
