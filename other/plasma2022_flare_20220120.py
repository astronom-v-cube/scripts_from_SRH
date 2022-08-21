#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 03:09:40 2022

@author: sergey_lesovoi
"""

import pylab as PL
from datetime import datetime
from astropy.io import fits
from optparse import OptionParser
import os
import urllib.request
from scipy.io import readsav
import base2uvw_612
from astropy import constants
import numpy as NP

def hm_format(t, pos):
  global srhT0
  global srhdT
  T = t*srhdT + srhT0    
  hh = int(T / 3600.)
  T -= hh*3600.
  mm = int(T / 60.)
  return '%02d:%02d' % (hh,mm);

def freq_format(f, pos):
  global srhFreqList
  if (f >= 0 and f < srhFreqList.shape[0]):
      return '%.1f' % (srhFreqList[int(f)]*1e-6)
  else:
      return '0'

def uv_format(uv, pos):
    global UVmax
    global UV
    return '%d' % (uv*(UVmax/UV))

def size_format(s, pos):
    return '%.1f' % (s*100)

srhCpFitPath = 'srh_0612_cp_20220120.fits'

srhCpFitFile = fits.open(srhCpFitPath)

t0 = 4200
t1 = 5500
dateList = srhCpFitFile[0].header['DATE-OBS'].split('-')
srhFreqList = srhCpFitFile[1].data['frequencies']
srhTime = srhCpFitFile[2].data['time']
srhT0 = srhTime[0,t0]
srhdT = srhTime[0,1] - srhTime[0,0]
srhCorrI = srhCpFitFile[2].data['I']
srhCorrV = srhCpFitFile[2].data['V']
srhFluxI = srhCpFitFile[2].data['flux_I']
srhFluxV = srhCpFitFile[2].data['flux_V']

fig = PL.figure(figsize=(12,6))
pl = fig.subplots(nrows=2,ncols=2)
for rr in range(2):
    for cc in range(2):
        pl[rr,cc].xaxis.set_major_formatter(PL.FuncFormatter(hm_format))
        #pl[0,0].xaxis.set_major_locator(PL.MultipleLocator(120))
        pl[rr,cc].yaxis.set_major_formatter(PL.FuncFormatter(freq_format))

pos0 = pl[0,0].imshow(srhFluxI[:,t0:t1],vmax=800,aspect=50,cmap='gnuplot')
pos1 = pl[1,0].imshow(srhCorrI[:,t0:t1]*100,vmax=30,aspect=50,cmap='gnuplot')
pos2 = pl[0,1].imshow(srhFluxV[:,t0:t1],vmax=150,aspect=50,cmap='gnuplot')
pos3 = pl[1,1].imshow(srhCorrV[:,t0:t1]*100,vmax=5,aspect=50,cmap='gnuplot')

pl[0,0].set_xlabel('UT')
pl[0,0].set_ylabel('GHz')
cb0 = fig.colorbar(pos0,ax=pl[0,0])
cb0.ax.set_title('s.f.u.')

pl[1,0].set_xlabel('UT')
pl[1,0].set_ylabel('GHz')
cb1 = fig.colorbar(pos1,ax=pl[1,0])
cb1.ax.set_title('%')

pl[0,1].set_xlabel('UT')
pl[0,1].set_ylabel('GHz')
cb2 = fig.colorbar(pos2,ax=pl[0,1])
cb2.ax.set_title('s.f.u.')

pl[1,1].set_xlabel('UT')
pl[1,1].set_ylabel('GHz')
cb3 = fig.colorbar(pos3,ax=pl[1,1])
cb3.ax.set_title('%')

startTime = phaseEdit.srhFits.freqTime[0,0]
timeInd = NP.where(srhTime[0] == startTime)[0][0]
T = phaseEdit.srhFits.visLcp.shape[1]
UV = 1500
UVmax = 15e3
deltaUV = UVmax / UV
dynUvSpectrum = NP.zeros((T,UV))

lambdas = (constants.c / (phaseEdit.srhFits.freqList * 1e3)).to_value()
delta = phaseEdit.srhFits.getDeclination()
for tt in range(T):
    hA = phaseEdit.srhFits.getHourAngle(tt)
    uvSpectrumSum = NP.ones(UV)
    for vv in range(8192):
#    for vv in range(phaseEdit.srhFits.visLcp.shape[2]):
        uvw = base2uvw_612.base2uvw(hA,delta,phaseEdit.srhFits.antennaA[vv] + 1,phaseEdit.srhFits.antennaB[vv] + 1)
        uvDist = NP.sqrt(uvw[0]**2 + uvw[1]**2)
        for ff in range(phaseEdit.srhFits.visLcp.shape[0]):
            uvDistInd = int((uvDist / lambdas[ff]) / deltaUV + .5)
            if (uvDistInd < UV):
                dynUvSpectrum[tt,uvDistInd] += (NP.abs(phaseEdit.srhFits.visLcp[ff,tt,vv]) + NP.abs(phaseEdit.srhFits.visRcp[ff,tt,vv]))
                uvSpectrumSum[uvDistInd] += 1
    dynUvSpectrum[tt] /= uvSpectrumSum
    print(tt)
        
dynUvSpectrumWidth = NP.zeros(T)
for tt in range(T):
    curMax = dynUvSpectrum[tt].max()
    try:
        halfInd = NP.where(NP.abs(curMax/2 - dynUvSpectrum[tt]) < curMax/50)
        dynUvSpectrumWidth[tt] = NP.mean(halfInd)
    except:
        dynUvSpectrumWidth[tt] = 1.

cmap='hot'
fig = PL.figure(figsize=(12,8))
pl = fig.subplots(nrows=3,ncols=1)

pl[0].xaxis.set_major_formatter(PL.FuncFormatter(hm_format))
pl[1].xaxis.set_major_formatter(PL.FuncFormatter(hm_format))
pl[2].xaxis.set_major_formatter(PL.FuncFormatter(hm_format))
pl[0].yaxis.set_major_formatter(PL.FuncFormatter(freq_format))
pl[1].yaxis.set_major_formatter(PL.FuncFormatter(freq_format))
pl[2].yaxis.set_major_formatter(PL.FuncFormatter(uv_format))

pos0 = pl[0].imshow(srhFluxI[:,timeInd:timeInd+460],vmax=900,aspect=7,cmap=cmap,origin='lower')
pos1 = pl[1].imshow(srhCorrI[:,timeInd:timeInd+460]*100,vmax=30,aspect=7,cmap=cmap,origin='lower')
pos2 = pl[2].imshow(dynUvSpectrum.T,vmax=.8,aspect=.075,cmap=cmap,origin='lower')

pl[0].set_xlabel('UT')
pl[0].set_ylabel('GHz')
cb0 = fig.colorbar(pos0,ax=pl[0])
cb0.ax.set_title('s.f.u.')

pl[1].set_xlabel('UT')
pl[1].set_ylabel('GHz')
cb1 = fig.colorbar(pos1,ax=pl[1])
cb1.ax.set_title('%')

pl[2].set_xlabel('UT')
pl[2].set_ylabel('uv_distance')
cb1 = fig.colorbar(pos1,ax=pl[2])
cb1.ax.set_title('au')

tpl = pl[0].twinx()
tpl.yaxis.set_major_formatter(PL.FuncFormatter(size_format))
tpl.plot(1./dynUvSpectrumWidth,'-',markersize=1.,color='green')
tpl.set_ylim(0,0.006)

fig=PL.figure(figsize=(6,4))
pl = fig.subplots(nrows=1,ncols=1)
pl.xaxis.set_major_formatter(PL.FuncFormatter(uv_format))
pl.plot(dynUvSpectrum[0:40].mean(axis=0),color='black',label='Solar disk response')
pl.set_xlim(0,100)
pl.set_xlabel('uv distance')
pl.grid()
pl.legend()
