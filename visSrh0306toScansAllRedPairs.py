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
import scipy.optimize as SPO
import scipy.signal as SPS

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

northPhaMatrix = MG.phaMatrixGen(31)
westPhaMatrix = MG.phaMatrixGen(32)
eastWestPhaMatrix = MG.phaMatrixGen(65)
ewnPhaMatrix = MG.phaMatrixGenEWN(65,31)

#path = '/home/sergeyvlesovoi/SRH36/20210206'
#path = 'SRH36_temp_20210206_2'
#path = 'SRH36_temp_20210222_2'
#path = 'SRH36_temp_20210223_4'
#path = 'SRH36_temp_20210224_1'
#path = 'SRH36_temp_20210224_3'
#path = 'SRH36_temp_20210225_2'
path = 'SRH36_temp_20210322_2'

fits_list = findFits(path,'*.fit')
fits_list.sort()

nfits = fits.open(fits_list[0])
time = nfits[1].data['time']
frequencyList = nfits[1].data['frequency']
samplesNumber = time.shape[1]
frequencyNumber = frequencyList.shape[0]
antennasNumber =  nfits[2].data['ant_name'].shape[0]
amplcp1 = nfits[1].data['amp_lcp']
amprcp1 = nfits[1].data['amp_rcp']
vislcp1 = nfits[1].data['vis_lcp']
visrcp1 = nfits[1].data['vis_rcp']
amplcp = amplcp1.reshape(frequencyNumber,samplesNumber,antennasNumber)
amprcp = amprcp1.reshape(frequencyNumber,samplesNumber,antennasNumber)
vislcp = vislcp1.reshape(frequencyNumber,samplesNumber,8128)
visrcp = visrcp1.reshape(frequencyNumber,samplesNumber,8128)

for file in fits_list[1:]:
    nfits = fits.open(file)
    time = NP.concatenate((time, nfits[1].data['time']), axis = 1)
    amplcp1 = nfits[1].data['amp_lcp']
    amprcp1 = nfits[1].data['amp_rcp']
    amplcp = NP.concatenate((amplcp,amplcp1.reshape(frequencyNumber,samplesNumber,antennasNumber)),axis=1)
    amprcp = NP.concatenate((amprcp,amprcp1.reshape(frequencyNumber,samplesNumber,antennasNumber)),axis=1)
    vislcp1 = nfits[1].data['vis_lcp']
    visrcp1 = nfits[1].data['vis_rcp']
    vislcp = NP.concatenate((vislcp, vislcp1.reshape(frequencyNumber,samplesNumber,8128)),axis=1)
    visrcp = NP.concatenate((visrcp, visrcp1.reshape(frequencyNumber,samplesNumber,8128)),axis=1)

ampScale = 1/(2e6*49*128)
PL.figure(figsize=(10,8))
pl0 = PL.subplot(111)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl0.set_title('SRH36 20210322, 2800, 4300 MHz')
pl0.plot(time[0,:],NP.abs(vislcp[3,:,0])*ampScale)

pl0.plot(time[0,:],NP.abs(vislcp[1,:,3009])*ampScale,label='N1 N2, 3100')
pl0.plot(time[0,:],NP.abs(vislcp[1,:,3008])*ampScale,label='N1 N2, 3100')
pl0.plot(time[0,:],NP.abs(vislcp[1,:,3007])*ampScale,label='N1 N2, 3100')

pl0.plot(time[0,:],NP.abs(vislcp[3,:,3009])*ampScale,label='N1 N2, 4300')
pl0.plot(time[0,:],NP.abs(vislcp[3,:,3008])*ampScale,label='N2 N3, 4300')
pl0.plot(time[0,:],NP.abs(vislcp[3,:,3007])*ampScale,label='N3 N4, 4300')

pl0.plot(time[0,:],NP.abs(vislcp[5,:,3009])*ampScale,label='N1 N2, 5600')
pl0.plot(time[0,:],NP.abs(vislcp[5,:,3008])*ampScale,label='N1 N2, 5600')
pl0.plot(time[0,:],NP.abs(vislcp[5,:,3007])*ampScale,label='N1 N2, 5600')
pl0.legend()
pl0.grid()

PL.figure(figsize=(10,8))
pl0 = PL.subplot(111)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl0.set_title('SRH36 20210322')
pl0.plot(time[0,:],vislcp[1,:,3008].real*ampScale,label='N1 N2 real, 3100')
pl0.plot(time[0,:],vislcp[1,:,3008].imag*ampScale,label='N1 N2 imag, 3100')
pl0.plot(time[0,:],NP.abs(vislcp[1,:,3008])*ampScale,label='N1 N2 amp, 3100')
pl0.plot(time[0,:],NP.angle(vislcp[1,:,3008])*0.0005,label='N1 N2 phase*0.0005, 3100')
pl0.legend()
pl0.grid()

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
#for scan in NP.linspace(0,samplesNumber-1,samplesNumber,dtype='int'):
#    for pair in westIndList:
#        vis_ftv_lcp[frequency,scan,pair] = NP.conj(vis_ftv_lcp[frequency,scan,pair])
#        vis_ftv_rcp[frequency,scan,pair] = NP.conj(vis_ftv_rcp[frequency,scan,pair])

northRedVis = NP.zeros(31,dtype='complex')
eastWestPhase = (NP.angle(vis_ftv_lcp[frequency,0,westInd0:westInd0+64]))
northRedVis[0]  = vis_ftv_lcp[frequency,0,32]
northRedVis[1:] = vis_ftv_lcp[frequency,0,northInd0:northInd0+30]
northPhase = (NP.angle(northRedVis))
ewnPhase = NP.concatenate((eastWestPhase,northPhase))
ewnAntPhase, c, d, e = NP.linalg.lstsq(ewnPhaMatrix,ewnPhase)
ewnAntPhaLcp = ewnAntPhase[2:]
northAntPhaLcp = ewnAntPhaLcp[65:]
eastWestAntPhaLcp = ewnAntPhaLcp[:65]

eastWestPhase = (NP.angle(vis_ftv_rcp[frequency,0,westInd0:westInd0+64]))
northRedVis[0]  = vis_ftv_rcp[frequency,0,32]
northRedVis[1:] = vis_ftv_rcp[frequency,0,northInd0:northInd0+30]
northPhase = (NP.angle(northRedVis))
ewnPhase = NP.concatenate((eastWestPhase,northPhase))
ewnAntPhase, c, d, e = NP.linalg.lstsq(ewnPhaMatrix,ewnPhase)
ewnAntPhaRcp = ewnAntPhase[2:]
northAntPhaRcp = ewnAntPhaRcp[65:]
eastWestAntPhaRcp = ewnAntPhaRcp[:65]

#northAntPhaLcpObj = NP.zeros(31)
#opt = SPO.minimize(objVisSum,northAntPhaLcpObj,complex_to_real(vis_ftv_lcp[frequency,0,northInd0:northInd0+30]))
#northAntPhaLcp = opt.x
#northAntPhaRcpObj = NP.zeros(31)
#opt = SPO.minimize(objVisSum,northAntPhaRcpObj,complex_to_real(vis_ftv_rcp[frequency,0,northInd0:northInd0+30]))
#northAntPhaRcp = opt.x
#
#westAntPhaLcpObj = NP.zeros(32)
#opt = SPO.minimize(objVisSum,westAntPhaLcpObj,complex_to_real(vis_ftv_lcp[frequency,0,westIndList[0:31]]))
#westAntPhaLcp = opt.x
#westAntPhaRcpObj = NP.zeros(32)
#opt = SPO.minimize(objVisSum,westAntPhaRcpObj,complex_to_real(vis_ftv_rcp[frequency,0,westIndList[0:31]]))
#westAntPhaRcp = opt.x
#
#eastAntPhaLcpObj = NP.zeros(32)
#opt = SPO.minimize(objVisSum,eastAntPhaLcpObj,complex_to_real(vis_ftv_lcp[frequency,0,eastIndList[0:31]]))
#eastAntPhaLcp = opt.x
#eastAntPhaRcpObj = NP.zeros(32)
#opt = SPO.minimize(objVisSum,eastAntPhaRcpObj,complex_to_real(vis_ftv_rcp[frequency,0,eastIndList[0:31]]))
#eastAntPhaRcp = opt.x

#eastWestAntPhaLcpObj = NP.zeros(65)
#opt = SPO.minimize(objVisSum,eastWestAntPhaLcpObj,complex_to_real(vis_ftv_lcp[frequency,0,westInd0:westInd0+64]))
#eastWestAntPhaLcp = opt.x
#eastWestAntPhaRcpObj = NP.zeros(65)
#opt = SPO.minimize(objVisSum,eastWestAntPhaRcpObj,complex_to_real(vis_ftv_rcp[frequency,0,westInd0:westInd0+64]))
#eastWestAntPhaRcp = opt.x

northScansLcp = []
northScansRcp = []
westScansRcp = []
westScansLcp = []
eastScansRcp = []
eastScansLcp = []
eastWestScansLcp = []
eastWestScansRcp = []

northAntPhaLcp = SPS.detrend(northAntPhaLcp)
northAntPhaRcp = SPS.detrend(northAntPhaRcp)
westAntPhaLcp = SPS.detrend(westAntPhaLcp)
westAntPhaRcp = SPS.detrend(westAntPhaRcp)
eastAntPhaLcp = SPS.detrend(eastAntPhaLcp)
eastAntPhaRcp = SPS.detrend(eastAntPhaRcp)
#eastWestAntPhaLcp = SPS.detrend(eastWestAntPhaLcp)
#eastWestAntPhaRcp = SPS.detrend(eastWestAntPhaRcp)

#northAntPhaLcpObj1D = northAntPhaLcp
#opt = SPO.minimize(objVis1DScanSum,northAntPhaLcpObj1D,complex_to_real(vis_ftv_lcp[frequency,0,northInd0:northInd0+int((30+1)/2*30)]),method='Powell')
#northAntPhaLcp1D = opt.x
#northAntPhaRcpObj1D = northAntPhaRcp
#opt = SPO.minimize(objVis1DScanSum,northAntPhaRcpObj1D,complex_to_real(vis_ftv_rcp[frequency,0,northInd0:northInd0+int((30+1)/2*30)]),method='Powell')
#northAntPhaRcp1D = opt.x
#
#westAntPhaLcpObj1D = westAntPhaLcp
#opt = SPO.minimize(objVis1DScanSum,westAntPhaLcpObj1D,complex_to_real(vis_ftv_lcp[frequency,0,westIndList]),method='Powell')
#westAntPhaLcp1D = opt.x
#westAntPhaRcpObj1D = westAntPhaRcp
#opt = SPO.minimize(objVis1DScanSum,westAntPhaRcpObj1D,complex_to_real(vis_ftv_rcp[frequency,0,westIndList]),method='Powell')
#westAntPhaRcp1D = opt.x
#
#eastAntPhaLcpObj1D = eastAntPhaLcp
#opt = SPO.minimize(objVis1DScanSum,eastAntPhaLcpObj1D,complex_to_real(vis_ftv_lcp[frequency,0,eastIndList]),method='Powell')
#eastAntPhaLcp1D = opt.x
#eastAntPhaRcpObj1D = eastAntPhaRcp
#opt = SPO.minimize(objVis1DScanSum,eastAntPhaRcpObj1D,complex_to_real(vis_ftv_rcp[frequency,0,eastIndList]),method='Powell')
#eastAntPhaRcp1D = opt.x

#eastWestAntPhaLcpObj1D = NP.zeros(65)
#eastWestAntPhaLcpObj1D = eastWestAntPhaLcp
#opt = SPO.minimize(objVis1DScanSum,eastWestAntPhaLcpObj1D,complex_to_real(vis_ftv_lcp[frequency,0,eastWestIndList]),method='Powell')
#eastWestAntPhaLcp1D = opt.x
#eastWestAntPhaRcpObj1D = NP.zeros(65)
#eastWestAntPhaRcpObj1D = eastWestAntPhaRcp
#opt = SPO.minimize(objVis1DScanSum,eastWestAntPhaRcpObj1D,complex_to_real(vis_ftv_rcp[frequency,0,eastWestIndList]),method='Powell')
#eastWestAntPhaRcp1D = opt.x

#northAntPhaLcp1D = SPS.detrend(northAntPhaLcp1D)
#northAntPhaRcp1D = SPS.detrend(northAntPhaRcp1D)
#
#westAntPhaLcp1D = SPS.detrend(westAntPhaLcp1D)
#westhAntPhaRcp1D = SPS.detrend(westAntPhaRcp1D)
#
#eastAntPhaLcp1D = SPS.detrend(eastAntPhaLcp1D)
#easthAntPhaRcp1D = SPS.detrend(eastAntPhaRcp1D)

#eastWestAntPhaLcp1D = SPS.detrend(eastWestAntPhaLcp1D)
#eastWestAntPhaRcp1D = SPS.detrend(eastWestAntPhaRcp1D)

northAntPhaLcp1D = northAntPhaLcp
northAntPhaRcp1D = northAntPhaRcp

westAntPhaLcp1D = westAntPhaLcp
westAntPhaRcp1D = westAntPhaRcp

eastAntPhaLcp1D = eastAntPhaLcp
eastAntPhaRcp1D = eastAntPhaRcp

for scan in NP.linspace(0,samplesNumber-1,samplesNumber,dtype='int'):
    northScansLcp.append(build1DScan(vis_ftv_lcp[frequency,scan,northInd0:northInd0+int((30+1)/2*30)], 30, northAntPhaLcp1D).real)
    northScansRcp.append(build1DScan(vis_ftv_rcp[frequency,scan,northInd0:northInd0+int((30+1)/2*30)], 30, northAntPhaRcp1D).real)
    westScansLcp.append(build1DScan(vis_ftv_lcp[frequency,scan,westIndList], 31, westAntPhaLcp1D).real)
    westScansRcp.append(build1DScan(vis_ftv_rcp[frequency,scan,westIndList], 31, westAntPhaRcp1D).real)
    eastScansLcp.append(build1DScan(vis_ftv_lcp[frequency,scan,eastIndList], 31, eastAntPhaLcp1D).real)
    eastScansRcp.append(build1DScan(vis_ftv_rcp[frequency,scan,eastIndList], 31, eastAntPhaRcp1D).real)
#    eastWestScansLcp.append(build1DScan(vis_ftv_lcp[frequency,scan,eastWestIndList], 64, eastWestAntPhaLcp1D).real)
#    eastWestScansRcp.append(build1DScan(vis_ftv_rcp[frequency,scan,eastWestIndList], 64, eastWestAntPhaRcp1D).real)

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
#eastWestScansLcp = NP.array(eastWestScansLcp)
#eastWestScansRcp = NP.array(eastWestScansRcp)
#eastWestScansLcp -= eastWestScansLcp.min()
#eastWestScansRcp -= eastWestScansRcp.min()

#PL.figure(figsize=(10,8))
#pl0 = PL.subplot(111)
#pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
#pl0.set_title('SRH36 20210203, %d MHz'%(frequencyList[frequency]*1e-3))
#pl0.plot(northScansLcp.mean(0))
#pl0.plot(northScansRcp.mean(0))
#pl0.plot(westScansLcp.mean(0))
#pl0.plot(westScansRcp.mean(0))
#pl0.plot(eastScansLcp.mean(0))
#pl0.plot(eastScansRcp.mean(0))

PL.figure()
pl1 = PL.subplot(111)
pl1.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
pl1.set_title('SRH36 20210203, %d MHz'%(frequencyList[frequency]*1e-3))
pl1.plot(time[frequency],NP.abs(vislcp[frequency,:,westInd0 + 0]))
pl1.plot(time[frequency],NP.abs(vislcp[frequency,:,westInd0 + 1]))
pl1.plot(time[frequency],NP.abs(vislcp[frequency,:,westInd0 + 3]))
pl1.plot(time[frequency],NP.abs(vislcp[frequency,:,westInd0 + 4]))
pl1.set_ylim(0,2e7)
pl1.set_xlabel('UTC')

PL.figure()
pl2 = PL.subplot(111)
pl2.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
pl2.set_title('SRH36 20210203, 2800-3250 MHz')
for f in range(16):
    pl2.plot(time[f],NP.abs(vislcp[f,:,:]).sum(1))
pl2.set_ylim(0,5e9)
pl2.set_xlabel('UTC')


