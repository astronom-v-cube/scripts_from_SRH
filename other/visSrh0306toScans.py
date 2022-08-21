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
import phaMatrixGen as MG

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

def hhmm_format(t, pos):
  hh = int(t / 3600.)
  t -= hh*3600.
  mm = int(t / 60.)
  t -= mm*60
  ss = int(t)
  return '%02d:%02d:%02d' % (hh,mm,ss)

northPhaMatrix = MG.phaMatrixGen(31)
westPhaMatrix = MG.phaMatrixGen(32)
eastWestPhaMatrix = MG.phaMatrixGen(65)
ewnPhaMatrix = MG.phaMatrixGenEWN(65,31)

path = 'SRH36_temp_20210320_1'
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
eastInd0 = 3505
westIndList = []
eastIndList = []
eastWestIndList = []
eastWestInd = 0
westInd = 0
northInd0 = 3007

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

northVsrc = NP.zeros(512,dtype='complex')
northV_ = NP.zeros(512,dtype='complex')

vis_ftv_lcp = vislcp
vis_ftv_rcp = visrcp
frequency = 3
for scan in NP.linspace(0,samplesNumber-1,samplesNumber,dtype='int'):
    for pair in westIndList:
        vis_ftv_lcp[frequency,scan,pair] = NP.conj(vis_ftv_lcp[frequency,scan,pair])
        vis_ftv_rcp[frequency,scan,pair] = NP.conj(vis_ftv_rcp[frequency,scan,pair])

northRedVis = NP.zeros(31,dtype='complex')
#eastWestPhase = NP.unwrap(NP.angle(vis_ftv_lcp[frequency,0,westInd0:westInd0+64]))
eastWestPhase = (NP.angle(vis_ftv_lcp[frequency,0,westInd0:westInd0+64]))
northRedVis[0]  = vis_ftv_lcp[frequency,0,32]
northRedVis[1:] = vis_ftv_lcp[frequency,0,northInd0:northInd0+30]
#northPhase = NP.unwrap(NP.angle(northRedVis))
northPhase = (NP.angle(northRedVis))
ewnPhase = NP.concatenate((eastWestPhase,northPhase))
ewnAntPhase, c, d, e = NP.linalg.lstsq(ewnPhaMatrix,ewnPhase)
ewnAntPhaLcp = ewnAntPhase[2:]
northAntPhaLcp = ewnAntPhaLcp[65:]
westAntPhaLcp = ewnAntPhaLcp[:32]
eastAntPhaLcp = ewnAntPhaLcp[33:65]
eastWestAntPhaLcp = ewnAntPhaLcp[:65]

#eastWestPhase = NP.unwrap(NP.angle(vis_ftv_rcp[frequency,0,westInd0:westInd0+64]))
eastWestPhase = (NP.angle(vis_ftv_rcp[frequency,0,westInd0:westInd0+64]))
northRedVis[0]  = vis_ftv_rcp[frequency,0,32]
northRedVis[1:] = vis_ftv_rcp[frequency,0,northInd0:northInd0+30]
#northPhase = NP.unwrap(NP.angle(northRedVis))
northPhase = (NP.angle(northRedVis))
ewnPhase = NP.concatenate((eastWestPhase,northPhase))
ewnAntPhase, c, d, e = NP.linalg.lstsq(ewnPhaMatrix,ewnPhase)
ewnAntPhaRcp = ewnAntPhase[2:]
northAntPhaRcp = ewnAntPhaRcp[65:]
westAntPhaRcp = ewnAntPhaRcp[:32]
eastAntPhaRcp = ewnAntPhaRcp[33:65]
eastWestAntPhaRcp = ewnAntPhaRcp[:65]

northScansLcp = []
northScansRcp = []
westScansLcp = []
westScansRcp = []
eastScansLcp = []
eastScansRcp = []
eastWestScansLcp = []
eastWestScansRcp = []

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

    eastInd = 0
    eastU_ = NP.zeros(512,dtype='complex')
    eastVis_lcp = vis_ftv_lcp[frequency,scan,eastIndList]
    for i in range(31):
        for j in range(31 - i):
            eastU_[i + 1] += eastVis_lcp[eastInd] * NP.exp(1j*(-eastAntPhaLcp[j] + eastAntPhaLcp[j+i+1]))
            eastInd += 1
        eastU_[512 - i - 1] = NP.conj(eastU_[i+1])
        
    eastScan = NP.roll(NP.fft.fft(eastU_),256)
    eastScansLcp.append(eastScan.real)

    eastInd = 0
    eastU_ = NP.zeros(512,dtype='complex')
    eastVis_rcp = vis_ftv_rcp[frequency,scan,eastIndList]
    for i in range(31):
        for j in range(31 - i):
            eastU_[i + 1] += eastVis_rcp[eastInd] * NP.exp(1j*(-eastAntPhaLcp[j] + eastAntPhaLcp[j+i+1]))
            eastInd += 1
        eastU_[512 - i - 1] = NP.conj(eastU_[i+1])
        
    eastScan = NP.roll(NP.fft.fft(eastU_),256)
    eastScansRcp.append(eastScan.real)

    eastWestInd = 0
    eastWestSpectrum = NP.zeros(512,dtype='complex')
    eastWestVis_lcp = vis_ftv_lcp[frequency,scan,eastWestIndList]
    for i in range(64):
        for j in range(64 - i):
            eastWestSpectrum[i + 1] += eastWestVis_lcp[eastWestInd] * NP.exp(1j*(-eastWestAntPhaLcp[j] + eastWestAntPhaLcp[j+i+1]))
            eastWestInd += 1
        eastWestSpectrum[512 - i - 1] = NP.conj(eastWestSpectrum[i+1])
        
    eastWestScan = NP.roll(NP.fft.fft(eastWestSpectrum),256)
    eastWestScansLcp.append(eastWestScan.real)

    eastWestInd = 0
    eastWestSpectrum = NP.zeros(512,dtype='complex')
    eastWestVis_rcp = vis_ftv_rcp[frequency,scan,eastWestIndList]
    for i in range(64):
        for j in range(64 - i):
            eastWestSpectrum[i + 1] += eastWestVis_rcp[eastWestInd] * NP.exp(1j*(-eastWestAntPhaRcp[j] + eastWestAntPhaRcp[j+i+1]))
            eastWestInd += 1
        eastWestSpectrum[512 - i - 1] = NP.conj(eastWestSpectrum[i+1])
        
    eastWestScan = NP.roll(NP.fft.fft(eastWestSpectrum),256)
    eastWestScansRcp.append(eastWestScan.real)

northScansLcp = NP.array(northScansLcp)
northScansRcp = NP.array(northScansRcp)
northScansLcp -= northScansLcp.min()
northScansRcp -= northScansRcp.min()
westScansLcp = NP.array(westScansLcp)
westScansRcp = NP.array(westScansRcp)
westScansLcp -= westScansLcp.min()
westScansRcp -= westScansRcp.min()
eastScansLcp = NP.array(eastScansLcp)
eastScansRcp = NP.array(eastScansRcp)
eastScansLcp -= eastScansLcp.min()
eastScansRcp -= eastScansRcp.min()
eastWestScansLcp = NP.array(eastWestScansLcp)
eastWestScansRcp = NP.array(eastWestScansRcp)
eastWestScansLcp -= eastWestScansLcp.min()
eastWestScansRcp -= eastWestScansRcp.min()
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
PL.figure(figsize=(10,8))
pl0 = PL.subplot(111)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
#pl0.set_title('SRH36 20210203, 3 GHz')
pl0.set_xlabel('UT')
pl0.set_ylabel('arbitrary')
#pl0.plot(time[frequency,:],northScansLcp[:,200],label='North')
#pl0.plot(time[frequency,:],eastWestScansLcp[:,200], label='EastWest')
pl0.plot(time[frequency,:],northScansRcp[:,200],label='North')
pl0.plot(time[frequency,:],eastWestScansRcp[:,200], label='EastWest')
#pl0.plot(time[frequency,:],westScansLcp[:,200] + westScansRcp[:,200],label='West')
#pl0.plot(time[frequency,:],eastScansLcp[:,200] + eastScansRcp[:,200],label='East')
pl0.legend()

PL.figure(figsize=(10,8))
pl1 = PL.subplot(111)
pl1.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl1.plot(time[frequency,:], NP.abs(vislcp[frequency,:,:]).sum(1))
pl1.plot(time[frequency,:], NP.abs(visrcp[frequency,:,:]).sum(1))

