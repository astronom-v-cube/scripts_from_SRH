#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 08:22:44 2021

@author: svlesovoi
"""
import numpy as NP
import pylab as PL
from astropy.io import fits
import os
import fnmatch

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

northPhaMatrix = [ \
            [1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1]]

westPhaMatrix = [ \
            [1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1]]

#rawData = fits.open('srh_20210108T030001.fit')
#rawData = fits.open('srh_20210108T021800.fit')
#rawData = fits.open('srh_20210108T044236.fit')
#rawData = fits.open('srh_20210108T044301.fit')
#rawData = fits.open('srh_20210108T060730.fit')
#rawData = fits.open('srh_20210108T061118.fit')
#rawData = fits.open('srh_20210108T061355.fit')
#rawData = fits.open('srh_20210108T061642.fit')
#rawData = fits.open('srh_20210108T061843.fit')
#rawData = fits.open('srh_20210108T062137.fit')
#rawData = fits.open('srh_20210108T062345.fit')
#rawData = fits.open('srh_20210108T063718.fit')
#rawData = fits.open('srh_20210108T064345.fit')
#rawData = fits.open('srh_20210108T070505.fit')
#rawData = fits.open('srh_20210108T071103.fit')
#rawData = fits.open('srh_20210108T071129.fit')
#rawData = fits.open('srh_20210108T071154.fit')
#rawData = fits.open('srh_20210108T071220.fit')
#rawData = fits.open('srh_20210108T071245.fit')
#rawData = fits.open('srh_20210108T071311.fit')
#rawData = fits.open('srh_20210108T071337.fit')

#path = 'SRH36_temp'
#fits_list = findFits(path,'*20210108T071*.fit')
#path = 'SRH36_temp_20210109'
#fits_list = findFits(path,'*.fit')
#path = 'SRH36_temp_20210109_1'
#path = 'SRH36_temp_20210109_2'
#path = 'SRH36_temp_20210109_3'
#path = 'SRH36_temp_20210202_1'
path = 'SRH36_temp_20210203_1'
fits_list = findFits(path,'*.fit')
fits_list.sort()

nfits = fits.open(fits_list[0])
time = nfits[1].data['time']
freqList = nfits[1].data['frequency']
samplesNumber = time.shape[1]
freqeuncyNumber = freqList.shape[0]
vislcp1 = nfits[1].data['vis_lcp']
visrcp1 = nfits[1].data['vis_rcp']
vislcp = vislcp1.reshape(freqeuncyNumber,samplesNumber,8128)
visrcp = visrcp1.reshape(freqeuncyNumber,samplesNumber,8128)

for file in fits_list[1:]:
    nfits = fits.open(file)
    time = NP.concatenate((time, nfits[1].data['time']), axis = 1)
    vislcp1 = nfits[1].data['vis_lcp']
    visrcp1 = nfits[1].data['vis_rcp']
    vislcp = NP.concatenate((vislcp, vislcp1.reshape(freqeuncyNumber,samplesNumber,8128)),axis=1)
    visrcp = NP.concatenate((visrcp, visrcp1.reshape(freqeuncyNumber,samplesNumber,8128)),axis=1)

samplesNumber = time.shape[1]
westInd0 = 3472
westIndList = []
eastWestInd = 0
westInd = 0
northInd0 = 3007

for i in range(31):
    for j in range(31 - i):
        westIndList.append(westInd0 + j + eastWestInd)
    eastWestInd += 96 - i

northVsrc = NP.zeros(512,dtype='complex')
northV_ = NP.zeros(512,dtype='complex')

vis_ftv_lcp = vislcp
vis_ftv_rcp = visrcp
frequency = 2

northPhase = NP.unwrap(NP.angle(vis_ftv_lcp[frequency,0,northInd0:northInd0+30]))
northAntPhaLcp, c, d, e = NP.linalg.lstsq(northPhaMatrix,northPhase)
northAntPhaLcp = northAntPhaLcp[1:]

northPhase = NP.unwrap(NP.angle(vis_ftv_rcp[frequency,0,northInd0:northInd0+30]))
northAntPhaRcp, c, d, e = NP.linalg.lstsq(northPhaMatrix,northPhase)
northAntPhaRcp = northAntPhaRcp[1:]

westPhase = NP.unwrap(NP.angle(vis_ftv_lcp[frequency,0,westInd0:westInd0+31]))
westAntPhaLcp, c, d, e = NP.linalg.lstsq(westPhaMatrix,westPhase)
westAntPhaLcp = westAntPhaLcp[1:]

westPhase = NP.unwrap(NP.angle(vis_ftv_rcp[frequency,0,westInd0:westInd0+31]))
westAntPhaRcp, c, d, e = NP.linalg.lstsq(westPhaMatrix,westPhase)
westAntPhaRcp = westAntPhaRcp[1:]

northScansLcp = []
westScansLcp = []
northScansRcp = []
westScansRcp = []
for scan in NP.linspace(0,samplesNumber-1,samplesNumber,dtype='int'):
    northInd = 0
    northVsrc = NP.zeros(512,dtype='complex')
    northV_ = NP.zeros(512,dtype='complex')
    for i in range(30):
        for j in range(30 - i):
            northVsrc[i + 1] += vis_ftv_lcp[frequency,scan,northInd0 + northInd]
            northV_[i + 1] += vis_ftv_lcp[frequency,scan,northInd0 + northInd] * NP.exp(1j*(-northAntPhaLcp[j] + northAntPhaLcp[j+i+1]))
            northInd += 1
        northVsrc[512 - i - 1] = NP.conj(northVsrc[i+1])
        northV_[512 - i - 1] = NP.conj(northV_[i+1])
        
    northScan = NP.roll(NP.fft.fft(northV_),256)
    northScansLcp.append(northScan.real)
    
    northInd = 0
    northVsrc = NP.zeros(512,dtype='complex')
    northV_ = NP.zeros(512,dtype='complex')
    for i in range(30):
        for j in range(30 - i):
            northVsrc[i + 1] += vis_ftv_rcp[frequency,scan,northInd0 + northInd]
            northV_[i + 1] += vis_ftv_rcp[frequency,scan,northInd0 + northInd] * NP.exp(1j*(-northAntPhaRcp[j] + northAntPhaRcp[j+i+1]))
            northInd += 1
        northVsrc[512 - i - 1] = NP.conj(northVsrc[i+1])
        northV_[512 - i - 1] = NP.conj(northV_[i+1])
        
    northScan = NP.roll(NP.fft.fft(northV_),256)
    northScansRcp.append(northScan.real)

    westInd = 0
    westUsrc = NP.zeros(512,dtype='complex')
    westU_ = NP.zeros(512,dtype='complex')
    westVis_lcp = vis_ftv_lcp[frequency,scan,westIndList]
    for i in range(31):
        for j in range(31 - i):
            westUsrc[i + 1] += westVis_lcp[westInd]
            westU_[i + 1] += westVis_lcp[westInd] * NP.exp(1j*(-westAntPhaLcp[j] + westAntPhaLcp[j+i+1]))
            westInd += 1
        westUsrc[512 - i - 1] = NP.conj(westUsrc[i+1])
        westU_[512 - i - 1] = NP.conj(westU_[i+1])
        
    westScan = NP.roll(NP.fft.fft(westU_),256)
    westScansLcp.append(westScan.real)

    westInd = 0
    westUsrc = NP.zeros(512,dtype='complex')
    westU_ = NP.zeros(512,dtype='complex')
    westVis_rcp = vis_ftv_rcp[frequency,scan,westIndList]
    for i in range(31):
        for j in range(31 - i):
            westUsrc[i + 1] += westVis_lcp[westInd]
            westU_[i + 1] += westVis_rcp[westInd] * NP.exp(1j*(-westAntPhaLcp[j] + westAntPhaLcp[j+i+1]))
            westInd += 1
        westUsrc[512 - i - 1] = NP.conj(westUsrc[i+1])
        westU_[512 - i - 1] = NP.conj(westU_[i+1])
        
    westScan = NP.roll(NP.fft.fft(westU_),256)
    westScansRcp.append(westScan.real)

northScansLcp = NP.array(northScansLcp)
northScansRcp = NP.array(northScansRcp)
westScansLcp = NP.array(westScansLcp)
westScansRcp = NP.array(westScansRcp)
#PL.figure()
#PL.title('North LCP + RCP')
#PL.imshow(northScansLcp + northScansRcp)
#PL.figure()
#PL.title('North LCP - RCP')
#PL.imshow(northScansLcp - northScansRcp)
##
#PL.figure()
#PL.title('WestEast')
#PL.imshow(westScansLcp)
PL.figure()
PL.plot(time[frequency,:],northScansLcp[:,200] + northScansRcp[:,200])
PL.plot(time[frequency,:],westScansLcp[:,200] + westScansRcp[:,200])
