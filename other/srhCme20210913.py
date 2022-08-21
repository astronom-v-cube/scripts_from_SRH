#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 07:19:02 2021

@author: svlesovoi
"""

import numpy as NP
import pylab as PL;
from astropy.io import fits
import os, fnmatch, sys

def arcmin_format(xy, pos):
    global gArcsecPerPixel
    global gXSize
    return '%2d' % ((xy - gXSize/2) * gArcsecPerPixel / 60);

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

fitPath = 'cleanMaps/20210913/'
iFitNames =  findFits(fitPath,'*I*_2800.fit')
vFitNames =  findFits(fitPath,'*V*_2800.fit')
iFitNames.sort()
vFitNames.sort()
iImages = []
vImages = []
freqList = []
gXSize = 0
gArcsecPerPixel = 0
gDateObs = []
for fileName in iFitNames:
    fd = fits.open(fileName)
    iImages.append(NP.flipud(fd[0].data))
    freqList.append(fd[0].header['FREQUENC'])
    gXSize = int(fd[0].header['NAXIS1'])
    gArcsecPerPixel = fd[0].header['CDELT1']
    gDateObs.append(fd[0].header['DATE-OBS'].split('.')[0])
    fd.close()
      
for fileName in vFitNames:
    fd = fits.open(fileName)
    vImages.append(NP.flipud(fd[0].data))
    freqList.append(fd[0].header['FREQUENC'])
    gXSize = int(fd[0].header['NAXIS1'])
    gArcsecPerPixel = fd[0].header['CDELT1']
    gDateObs.append(fd[0].header['DATE-OBS'].split('.')[0])
    fd.close()

iImagesArr = NP.array(iImages)
dL = 100
trace = []
for l in range(dL):
    trace.append(iImagesArr[1:,90+l,340])
