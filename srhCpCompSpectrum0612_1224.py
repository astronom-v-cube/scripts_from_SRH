#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 03:25:49 2022

@author: sergey_lesovoi
"""

import numpy as NP
import pylab as PL
from matplotlib.ticker import MultipleLocator, FormatStrFormatter;
from astropy.io import fits
from astropy.time import Time, TimeDelta
import sunpy.image.resample as SUNIM;

def hm_format(t, pos):
    if (t < compLen):
        t = cpTime[int(t)]
        hh = int(t / 3600.)
        t -= hh*3600.
        mm = int(t / 60.)
        return '%02d:%02d' % (hh,mm);
    else:
        return '%02d:%02d' % (0,0);

def freq_format(f, pos):
    return '%2.1f' % (6 + f*0.4)

# fit0612 = fits.open('srh_0612_cp_20220314.fits')
# fit1224 = fits.open('srh_1224_cp_20220314.fits')
fit0612 = fits.open('srh_0612_cp_20220315.fits')
fit1224 = fits.open('srh_1224_cp_20220315.fits')

cpTime = fit0612[2].data['time'][0]

cpI0612 = fit0612[2].data['I']
cpI1224 = fit1224[2].data['I']
cpV0612 = fit0612[2].data['V']
cpV1224 = fit1224[2].data['V']

len0612 = cpI0612.shape[1]
len1224 = cpI1224.shape[1]

cpI1224 = SUNIM.resample(cpI1224,(cpI1224.shape[0]*2,len1224))
cpV1224 = SUNIM.resample(cpV1224,(cpV1224.shape[0]*2,len1224))

compLen = max([len0612,len1224])

cpComp = NP.zeros((48,compLen))
cpComp[0:16,0:len0612] = cpI0612
cpComp[16:48,0:len1224] = (cpI1224 - 0.0008)*7.5

fig = PL.figure()
fig.suptitle('SRH I 6-24 20220314')
pl = fig.subplots()
pl.xaxis.set_major_locator(MultipleLocator(600.));
pl.xaxis.set_major_formatter(PL.FuncFormatter(hm_format));
pl.yaxis.set_major_locator(MultipleLocator(8));
pl.yaxis.set_major_formatter(PL.FuncFormatter(freq_format));
pl.set_xlabel('UTC')
pl.set_ylabel('GHz')

pl.imshow(cpComp,aspect=100,cmap='jet',interpolation='bessel',vmin=0,vmax=0.03)

cpComp = NP.zeros((48,compLen))
cpComp[0:16,0:len0612] = cpV0612
cpComp[16:48,0:len1224] = cpV1224*7.5

fig = PL.figure()
fig.suptitle('SRH V 6-24 20220314')
pl = fig.subplots()
pl.xaxis.set_major_locator(MultipleLocator(600.));
pl.xaxis.set_major_formatter(PL.FuncFormatter(hm_format));
pl.yaxis.set_major_locator(MultipleLocator(8));
pl.yaxis.set_major_formatter(PL.FuncFormatter(freq_format));
pl.set_xlabel('UTC')
pl.set_ylabel('GHz')

pl.imshow(cpComp,aspect=100,cmap='gray_r',interpolation='bessel',vmin=-0.005,vmax=0.005)
