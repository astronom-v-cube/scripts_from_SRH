#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 16:29:25 2021

@author: sergeyvlesovoi
"""
from astropy.io import fits
import pylab as PL

loFits = fits.open('MF/20210331_mf_2800.fits')
hiFits = fits.open('MF/20210331_mf_5600.fits')
mfFits = fits.open('MF/20210331_mf_2800+5600.fits')

loData = loFits[0].data[0][0]
hiData = hiFits[0].data[0][0]
mfData = mfFits[0].data[0][0]

PL.figure()
PL.imshow(loData,origin='lower',cmap='hot')

PL.figure()
PL.imshow(hiData,origin='lower',cmap='hot')

PL.figure()
PL.imshow(mfData,origin='lower',cmap='hot')