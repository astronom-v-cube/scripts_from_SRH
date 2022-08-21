#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 15:45:09 2021

@author: mariagloba
"""

from astropy.io import fits
import numpy as NP
import pylab as PL
from astropy import coordinates
from astropy import constants
from BadaryRAO import BadaryRAO
from sunpy import coordinates as sunpy_coordinates
from skimage.transform import warp, AffineTransform
import matrices_conj as matrices
import phaMatrixGen as MG
import scipy.optimize as SPO
import scipy.signal as SPS
from scipy.stats import linregress
import os, fnmatch;

def antEWFormat(f, pos):
    if f < 64:
        return 'W%02d' % ((64 - f)//2)
    elif f == 64:
        return 'C00' % ((64 - f)//2)
    else:
        return 'E%02d' % ((f - 64)//2)

def antNFormat(f, pos):
    if f < 64:
        return 'N%02d' % ((64 - f)//2)
    elif f == 64:
        return 'C00' % ((64 - f)//2)
    else:
        return 'N%02d' % ((f - 64)//2)

def real_to_complex(z):
    return z[:len(z)//2] + 1j * z[len(z)//2:]

def complex_to_real(z):
    return NP.concatenate((NP.real(z), NP.imag(z)))

def objVisSum(phases, visibilities):
    z = real_to_complex(visibilities)
    res = 0 + 1j*0
    for k in range(len(z)):
        res += NP.exp(1j*(phases[k + 1] - phases[k])) * z[k]
    return -NP.abs(res)

def build1DScan(scanVis, visNumber, antPhases):
    scanSpectrum = NP.zeros(512,dtype='complex')
    pairInd = 0
    for pair in range(visNumber):
        for ant in range(visNumber - pair):
            scanSpectrum[pair + 1] += scanVis[pairInd] * NP.exp(1j*(antPhases[ant + 1 + pair] - antPhases[ant]))
            pairInd += 1
        scanSpectrum[512 - pair - 1] = NP.conj(scanSpectrum[pair + 1])
    return NP.roll(NP.fft.fft(scanSpectrum),256)
    
def objVis1DScanSum(phases, visibilities):
    return -NP.abs(build1DScan(real_to_complex(visibilities), len(phases) - 1, phases).real).sum()

def getPQScale(size, FOV, time, freq):
    cosP = NP.sin(hAngle) * NP.cos(RAO.declination)
    cosQ = NP.cos(hAngle) * NP.cos(RAO.declination) * NP.sin(RAO.observatory.lat) - NP.sin(RAO.declination) * NP.cos(RAO.observatory.lat)
    FOV_p = 2.*(constants.c / (freq*1e3)) / (RAO.base*NP.sqrt(1. - cosP**2.));
    FOV_q = 2.*(constants.c / (freq*1e3)) / (RAO.base*NP.sqrt(1. - cosQ**2.));
    
    return [int(size*FOV/FOV_p.to_value()), int(size*FOV/FOV_q.to_value())]
    
def getPQ2HDMatrix():
    gP =  NP.arctan(NP.tan(hAngle)*NP.sin(RAO.declination));
    gQ =  NP.arctan(-(NP.sin(RAO.declination) / NP.tan(hAngle) + NP.cos(RAO.declination) / (NP.sin(hAngle)*NP.tan(RAO.observatory.lat))));
    
    if hAngle > 0:
        gQ = NP.pi + gQ;
    g = gP - gQ;
      
    pqMatrix = NP.zeros((3,3))
    pqMatrix[0, 0] =  NP.cos(gP) - NP.cos(g)*NP.cos(gQ)
    pqMatrix[0, 1] = -NP.cos(g)*NP.cos(gP) + NP.cos(gQ)
    pqMatrix[1, 0] =  NP.sin(gP) - NP.cos(g)*NP.sin(gQ)
    pqMatrix[1, 1] = -NP.cos(g)*NP.sin(gP) + NP.sin(gQ)
    pqMatrix /= NP.sin(g)**2.
    pqMatrix[2, 2] = 1.
    return pqMatrix

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

fitNames = findFits('SRH36_temp_20210501_3','*.fit')

fitNames.sort()
fitsFile = fits.open(fitNames[-1])
#fitsFile = fits.open('SRH36_temp_20210318_1/srh_20210318T021631.fit')

visScale = 1/(2e6*49)
ampScale = visScale / 128
westInd0 = 3472
eastInd0 = 3505
northInd0 = 3007
westIndList = []
eastIndList = []
eastWestIndList = []
eastWestInd = 0
westInd = 0

for i in range(31):
    for j in range(31 - i):
        westIndList.append(westInd0 + j + eastWestInd)
    eastWestInd += 96 - i

eastWestInd = 0
for i in range(31):
    for j in range(31 - i):
        eastIndList.append(eastInd0 + j + eastWestInd)
    eastWestInd += 96 - i

eastWestInd = 0
for i in range(64):
    for j in range(64 - i):
        eastWestIndList.append(westInd0 + j + eastWestInd)
    eastWestInd += 96 - i

antZeroRow = fitsFile[3].data['ant_zero_row'][:97]
time = fitsFile[1].data['time']
date_obs = fitsFile[0].header['DATE-OBS']
time_obs = fitsFile[0].header['TIME-OBS']
freqListLength = fitsFile[1].data['frequency'].shape[0]
freqList = fitsFile[1].data['frequency']*1e3
samplesNumber = time.shape[1]
visibilityNumber = fitsFile[1].data['vis_lcp'].shape[1]//samplesNumber

visLcp = fitsFile[1].data['vis_lcp'].reshape(freqListLength,samplesNumber,visibilityNumber) * visScale
visRcp = fitsFile[1].data['vis_rcp'].reshape(freqListLength,samplesNumber,visibilityNumber) * visScale

scan = 10
freq = 5
polarization = 1

uvSize = 2049
uv = NP.zeros((uvSize, uvSize), dtype = 'complex')

if polarization == 0:
    for i in range(31):
        for j in range(65):
            O = uvSize//2
            uv[O + (i+1)*2, O + (j-32)*2] = visLcp[freq, scan, i*97+j]
            uv[O - (i+1)*2, O - (j-32)*2] = NP.conj(uv[O + (i+1)*2, O + (j-32)*2])
    for i in range(32):
        O = uvSize//2
        uv[O, O + (i-32)*2] = visLcp[freq, scan, antZeroRow[i]]
        uv[O, O + (32-i)*2] = NP.conj(uv[O, O + (i-32)*2])
else:
    for i in range(31):
        for j in range(65):
            O = uvSize//2
            uv[O + (i+1)*2, O + (j-32)*2] = visRcp[freq, scan, i*97+j]
            uv[O - (i+1)*2, O - (j-32)*2] = NP.conj(uv[O + (i+1)*2, O + (j-32)*2])
    for i in range(32):
        O = uvSize//2
        uv[O, O + (i-32)*2] = visRcp[freq, scan, antZeroRow[i]]
        uv[O, O + (32-i)*2] = NP.conj(uv[O, O + (i-32)*2])
        
image = NP.fft.fft2(NP.roll(NP.roll(uv,uvSize//2+1,0),uvSize//2+1,1));
image = NP.roll(NP.roll(image,uvSize//2-1,0),uvSize//2-1,1);
showUV = NP.abs(uv[uvSize//2-65:uvSize//2+65,uvSize//2-65:uvSize//2+65])
showUV = NP.clip(showUV,0,0.1*NP.max(showUV))
showUV = NP.roll(showUV,-1,axis=1)
showUV = NP.roll(showUV,-1,axis=0)

fig = PL.figure()
spUV = fig.add_subplot(1,1,1,title='UV %s, %s, %d MHz, %s'%(date_obs, time_obs, fitsFile[1].data['frequency'][freq]/1e3, 'RCP' if polarization  else 'LCP'))
spUV.set_xlabel('antenna EW')
spUV.set_ylabel('antenna N')
spUV.xaxis.set_major_locator(PL.MultipleLocator(16));
spUV.xaxis.set_minor_locator(PL.MultipleLocator(4));
spUV.xaxis.set_major_formatter(PL.FuncFormatter(antEWFormat))
spUV.yaxis.set_major_locator(PL.MultipleLocator(16));
spUV.yaxis.set_minor_locator(PL.MultipleLocator(4));
spUV.yaxis.set_major_formatter(PL.FuncFormatter(antNFormat))
spUV.imshow(showUV,cmap='hot')

