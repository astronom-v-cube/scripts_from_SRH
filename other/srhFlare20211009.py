#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 10:37:51 2021

@author: sergeyvlesovoi
"""
import netCDF4 as nc
import numpy as NP
import cftime
import pylab as PL
from datetime import datetime
from astropy.io import fits
from optparse import OptionParser
import os
import urllib.request
from scipy.io import readsav
import ZirinTb
import matplotlib


def hms_format(t, pos):
  hh = int(t / 3600.)
  t -= hh*3600.
  mm = int(t / 60.)
  t -= mm*60.
  ss = int(t)
  return '%02d:%02d:%02d' % (hh,mm,ss);

def hm_format(t, pos):
  global srhT0
  global srhdT
  T = t*srhdT + srhT0    
  hh = int(T / 3600.)
  T -= hh*3600.
  mm = int(T / 60.)
  return '%02d:%02d' % (hh,mm);

def freq_format(f, pos):
  global flareFreqs
  return '%.1f' % (flareFreqs*1e-6);

Zi = ZirinTb.ZirinTb()

srhCpFitFile = fits.open('srh_cp_20211009.fits')
flareFreqs = srhCpFitFile[1].data['frequencies']
flareTime = srhCpFitFile[2].data['time'][:,16800:20600]
flareFluxI = srhCpFitFile[2].data['flux_I'][:,16800:20600]
flareFluxV = srhCpFitFile[2].data['flux_V'][:,16800:20600]
flareCorrI = srhCpFitFile[2].data['I'][:,16800:20600]
flareCorrV = srhCpFitFile[2].data['V'][:,16800:20600]

srhT0 = flareTime[0,900]
srhdT = flareTime[0,1] - flareTime[0,0]
ind0 = 620
startFluxes = [73.5, 79.4, 93.7, 100, 100, 100]
for i in range(flareFreqs.shape[0]):
    startFluxes[i] = Zi.getSfuAtFrequncy(flareFreqs[i]*1e-6)/2.3
    
flareFluxICal = NP.zeros_like(flareFluxI)
flareFluxVCal = NP.zeros_like(flareFluxI)
for i in range(flareFluxI.shape[0]):
    flareFluxICal[i,:] = (flareFluxI[i,:] - flareFluxI[i,0]*.5)/flareFluxI[i,ind0]*startFluxes[i]
    flareFluxVCal[i,:] = (flareFluxV[i,:])/flareFluxI[i,ind0]*startFluxes[i]
for i in range(flareFluxI.shape[0]):
    flareFluxICal[i,:] -= flareFluxICal[i,0]

flareCorrICal = NP.zeros_like(flareCorrI)
flareCorrVCal = NP.zeros_like(flareCorrI)
for i in range(flareFluxI.shape[0]):
    flareCorrICal[i,:] = flareCorrI[i,:] - flareCorrI[i,0]
    flareCorrVCal[i,:] = flareCorrV[i,:]

c_list=matplotlib.colors.LinearSegmentedColormap.from_list(PL.cm.datad['gist_rainbow'], colors=['r','g','b'], N = 6)

fig = PL.figure(figsize=(8,6))
fig.suptitle('Flux ' + srhCpFitFile[0].header['DATE-OBS'])
pl = fig.subplots(nrows=1,ncols=1)
pl.xaxis.set_major_formatter(PL.FuncFormatter(hms_format))
pl.xaxis.set_major_locator(PL.MultipleLocator(300))

for i in range(flareFreqs.shape[0]):
    flareFluxICal[i,1598:] *= 1.1
    flareFluxVCal[i,1598:] *= 1.1
    pl.plot(flareTime[i],flareFluxICal[i],color=c_list(i),label='%.1f GHz'%(flareFreqs[i]*1e-6))
    pl.plot(flareTime[i],flareFluxVCal[i],color=c_list(i))
pl.set_ylabel("sfu")
pl.set_xlabel("Time [UT]")
pl.set_xlim(6.45*3600, 7.1*3600)
pl.legend()
pl.grid()

fig = PL.figure(figsize=(8,6))
fig.suptitle('Corrplot ' + srhCpFitFile[0].header['DATE-OBS'])
pl = fig.subplots(nrows=1,ncols=1)
pl.xaxis.set_major_formatter(PL.FuncFormatter(hms_format))
pl.xaxis.set_major_locator(PL.MultipleLocator(300))
for i in range(flareFreqs.shape[0]):
    flareCorrICal[i,1598:] *= 1.17
    flareCorrVCal[i,1598:] *= 1.17
    pl.plot(flareTime[i],flareCorrICal[i],color=c_list(i),label='%.1f GHz'%(flareFreqs[i]*1e-6))
    pl.plot(flareTime[i],flareCorrVCal[i],color=c_list(i))
pl.set_ylabel("correlation")
pl.set_xlabel("Time [UT]")
pl.set_xlim(6.45*3600, 7.1*3600)
pl.legend()
pl.grid()
