#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 01:51:42 2021

@author: svlesovoi
"""

from astropy.io import fits
import numpy as NP
import pylab as PL
from astropy import coordinates
from astropy import constants
from BadaryRAO import BadaryRAO
from scipy.optimize import least_squares
import sunpy.coordinates
import base2uvw_36
from skimage.transform import warp, AffineTransform
import scipy.signal
from zeep import Client
import os, fnmatch, sys

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

fitPath = '/var/run/media/svlesovoi/FTP/SRH/SRH0612/20210627/temp'
fitNames =  findFits(fitPath,'*.fit')
fitNames.sort()

sVisIndex0 = 8192
weVisIndex0 = sVisIndex0 + 2016

for fitName in fitNames:
    fd = fits.open(fitName)
    if (fitName == fitNames[0]):
        freqList = fd[1].data['frequency']
        time = fd[1].data['time']
        freqNumber = freqList.shape[0]
        scanNumber = time.shape[1]
        antennaNumber = fd[2].data['ant_name'].shape[0]
        visNumber = antennaNumber*(antennaNumber - 1) // 2
        ampLcp = fd[1].data['amp_lcp'].reshape(freqNumber,scanNumber,antennaNumber)
        ampRcp = fd[1].data['amp_rcp'].reshape(freqNumber,scanNumber,antennaNumber)
        visLcp = fd[1].data['vis_lcp'].reshape(freqNumber,scanNumber,visNumber)
        visRcp = fd[1].data['vis_rcp'].reshape(freqNumber,scanNumber,visNumber)
    else:
        ampLcp = NP.concatenate((ampLcp,fd[1].data['amp_lcp'].reshape(freqNumber,scanNumber,antennaNumber)),axis=1)
        ampRcp = NP.concatenate((ampRcp,fd[1].data['amp_rcp'].reshape(freqNumber,scanNumber,antennaNumber)),axis=1)
        visLcp = NP.concatenate((visLcp,fd[1].data['vis_lcp'].reshape(freqNumber,scanNumber,visNumber)),axis=1)
        visRcp = NP.concatenate((visRcp,fd[1].data['vis_rcp'].reshape(freqNumber,scanNumber,visNumber)),axis=1)

PL.figure()
pair = 4
#for freq in range(freqNumber):
for freq in range(1):
    PL.plot(visLcp[freq,:,weVisIndex0+pair].real,visLcp[freq,:,weVisIndex0+pair].imag,'.')
    PL.plot(visRcp[freq,:,weVisIndex0+pair].real,visRcp[freq,:,weVisIndex0+pair].imag,'.')

#PL.plot(NP.abs(visLcp[4,:,:]).sum(axis=0),'.')