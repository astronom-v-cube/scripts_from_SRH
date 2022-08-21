#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 02:44:52 2020

@author: svlesovoi
"""
from astropy.io import fits
import numpy as NP
import pylab as PL
import struct

nfits = fits.open('srh_20200620T010817.fit')
time = nfits[1].data['time']
samplesNumber = time.shape[1]
amplcp = nfits[1].data['amp_lcp']
amplitudeNumber = amplcp.shape[1]//samplesNumber
amplcp = amplcp.reshape(1,samplesNumber,amplitudeNumber)
rawlcp = nfits[1].data['vis_lcp']
visibilityNumber = rawlcp.shape[1]//samplesNumber
unpackFormat = ">%dl" % (rawlcp.shape[1]*2)
ff = NP.array(struct.unpack(unpackFormat,rawlcp))
fff = ff.reshape(rawlcp.shape[1],2)
fff_real = fff[:,0].reshape(1,samplesNumber,visibilityNumber)
fff_imag = fff[:,1].reshape(1,samplesNumber,visibilityNumber)

for ant in range(48):
    PL.plot(amplcp[0,:,ant])

#for pair in range(48):
#    PL.plot(fff_real[0,:,pair])
#    PL.plot(fff_imag[0,:,pair])