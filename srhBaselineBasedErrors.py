#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 14:11:02 2021

@author: sergeyvlesovoi
"""
import os, fnmatch;
#from srhFitsFile36_amp import SrhFitsFile
from srhFitsFile36 import SrhFitsFile
import numpy as NP;
import datetime as DT;
import pylab as PL;
import sys;
import matplotlib;
import phaMatrixGen

dt_major = 300.;
dt_minor = 30.;
freq_first = 2;
freq_intervals = 5;

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

dateName = DT.datetime.now().strftime("%Y%m%d");
#fitNames = findFits('/home/sergeyvlesovoi/SRH36/20210405','*.fit')
#fitNames = findFits('SRH36_temp_20210427_1','*.fit')
#fitNames = findFits('SRH36_temp_20210501_1','*.fit')
#fitNames = findFits('SRH36_temp_20210501_2','*.fit')
#fitNames = findFits('SRH36_temp_20210501_3','*.fit')
fitNames = findFits('SRH36_temp_20210504_1','*.fit')
fitNames.sort();
fitNames = fitNames[1:];

visScale = 1/(2e6*49)
ampScale = visScale / 128
firstBaseline = 3007
secondBaselineOffset = 30

for fName in fitNames:
    if (fName == fitNames[0]):
        sF = SrhFitsFile(fName, 257)
    else:
        sF.append(fName)

#fig = PL.figure()
#sub = fig.add_subplot(1,1,1);
#sub.set_ylabel('rad');
#sub.set_xlabel('UT');
#sub.xaxis.set_major_locator(PL.MultipleLocator(60));
#sub.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
#sub.xaxis.set_minor_locator(PL.MultipleLocator(10));
#sub.set_ylim(-3.5,3.5)
#for freq in range(sF.freqListLength):
#    PL.plot(sF.freqTime[freq], NP.unwrap(NP.angle(sF.visLcp[freq,:,closureBase])),'.',label='%d'%(sF.freqList[freq]*1e-3))
#sub.legend()

#freq = 1

#fig = PL.figure()
#sub = fig.add_subplot(1,1,1);
#sub.set_ylabel('rad');
#sub.set_xlabel('UT');
#sub.xaxis.set_major_locator(PL.MultipleLocator(60));
#sub.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
#sub.xaxis.set_minor_locator(PL.MultipleLocator(10));
#sub.set_ylim(-1,7)
#for freq in range(sF.freqListLength):
#    PL.plot(sF.freqTime[freq], NP.unwrap(NP.angle(sF.visLcp[freq,:,closureBase]) + NP.angle(sF.visLcp[freq,:,closureBase+1]) - NP.angle(sF.visLcp[freq,:,closureBase+30])),'.')
#    sub.plot(sF.freqTime[freq], (phaseA1A2),'.',label='A1A2')
#    sub.plot(sF.freqTime[freq], (phaseA2A3),'.',label='A2A3')
#    sub.plot(sF.freqTime[freq], (phaseA1A3),'.',label='A1A3')
#    sub.plot(sF.freqTime[freq], (phaseA1A2 + phaseA2A3 - phaseA1A3)%(2*NP.pi),'.',label='A1A2+A2A3-A1A3')
#sub.legend()

closurePhS = NP.zeros((sF.freqListLength,29))
closurePhEW = NP.zeros((sF.freqListLength,95))
clPhS = NP.zeros((sF.freqListLength,sF.freqTime.shape[1],29))
clPhEW = NP.zeros((sF.freqListLength,sF.freqTime.shape[1],95))
for freq in range(sF.freqListLength):
    firstBaseline = 3007
    secondBaselineOffset = 30
    for i in range(29):
        phaseA1A2 = NP.angle(sF.visLcp[freq,:,firstBaseline + i])
        phaseA2A3 = NP.angle(sF.visLcp[freq,:,firstBaseline + 1 + i])
        phaseA1A3 = NP.angle(sF.visLcp[freq,:,firstBaseline + secondBaselineOffset + i])
#        closurePhS[freq,i] = (((phaseA1A2 + phaseA2A3 - phaseA1A3)%(2*NP.pi)).mean())
        clPhS[freq,:,i] = phaseA1A2 + phaseA2A3 - phaseA1A3
        closurePhS[freq,i] = (((phaseA1A2 + phaseA2A3 - phaseA1A3)).mean())
    
    firstBaseline = 3472
    secondBaselineOffset = 96
    for i in range(95):
        phaseA1A2 = NP.angle(sF.visLcp[freq,:,firstBaseline + i])
        phaseA2A3 = NP.angle(sF.visLcp[freq,:,firstBaseline + 1 + i])
        phaseA1A3 = NP.angle(sF.visLcp[freq,:,firstBaseline + secondBaselineOffset + i])
#        closurePhEW[freq,i] = (((phaseA1A2 + phaseA2A3 - phaseA1A3)%(2*NP.pi)).mean())
        clPhEW[freq,:,i] = phaseA1A2 + phaseA2A3 - phaseA1A3
        closurePhEW[freq,i] = (((phaseA1A2 + phaseA2A3 - phaseA1A3)).mean())

#PL.figure()
#PL.title(sF.dateObs)
#for freq in range(sF.freqListLength):
##    PL.plot(NP.unwrap(closurePhS[freq]),'.',label='%d'%(sF.freqList[freq]*1e-3))
#    PL.plot((closurePhS[freq]),'.',label='%d'%(sF.freqList[freq]*1e-3))
#PL.legend()
#
#PL.figure()
#PL.title(sF.dateObs)
#for freq in range(sF.freqListLength):
##    PL.plot(NP.unwrap(closurePhEW[freq]),'.',label='%d'%(sF.freqList[freq]*1e-3))
#    PL.plot((closurePhEW[freq]),'.',label='%d'%(sF.freqList[freq]*1e-3))
#PL.legend()

phaseNoiseLevel = 0.1
nPhaseLcp = NP.angle(sF.visLcp[freq,:,3007:3007+30])
ewPhaseLcp = NP.angle(sF.visLcp[freq,:,3472:3472+96])
nPhaseNoiseLcp = NP.std(nPhaseLcp,axis=0)
ewPhaseNoiseLcp = NP.std(ewPhaseLcp,axis=0)

ewFlaggedVis = NP.array(NP.where(ewPhaseNoiseLcp > phaseNoiseLevel))
nFlaggedVis = NP.array(NP.where(nPhaseNoiseLcp > phaseNoiseLevel))
ewFlaggedVis = ewFlaggedVis[0][::2]
nFlaggedVis = nFlaggedVis[0][::2]

flaggedAnts = []
for k in range(ewFlaggedVis.shape[0]):
    antAName = str(sF.antennaA[ewFlaggedVis[k] + 3472])
    flaggedAnts.append(sF.hduList[2].data['ant_name'][NP.where(sF.antennaNumbers == antAName)])

for k in range(nFlaggedVis.shape[0]):
    antAName = str(sF.antennaA[nFlaggedVis[k] + 3007])
    flaggedAnts.append(sF.hduList[2].data['ant_name'][NP.where(sF.antennaNumbers == antAName)])

#sF.flagAntennasByPhaseNoise()

ewPhaMatrix = phaMatrixGen.phaMatrixGen(97)
antPhaEW = []
for i in range(sF.visLcp.shape[1]):
    visPha = ewPhaseLcp[i,:]
    #visPha[ewFlaggedVis] = 0.
    antPha, c, d, e = NP.linalg.lstsq(ewPhaMatrix,NP.angle(visPha))
    antPhaEW.append(antPha)

nPhaMatrix = phaMatrixGen.phaMatrixGen(31)
antPhaN = []
for i in range(sF.visLcp.shape[1]):
    visPha = nPhaseLcp[i,:]
#    visPha[nFlaggedVis] = 0.
    antPha, c, d, e = NP.linalg.lstsq(nPhaMatrix,NP.angle(visPha))
    antPhaN.append(antPha)
    
antPhaEW = NP.array(antPhaEW)
antPhaN = NP.array(antPhaN)

#PL.figure()
#PL.imshow(antPhaN,aspect = 0.2)

fig, pl = PL.subplots(nrows=3,ncols=1,figsize=(8,8))
fig.tight_layout()
pl[0].set_title('N1+N2-N3')
pl[0].text(sF.freqTime[0,1],360,'SRH closure phases %s'%(sF.dateObs.split('T')[0]),fontsize=14)
pl[1].set_title('N4+N5-N6')
pl[2].set_title('N7+N8-N9')
for freq in range(sF.freqListLength):
    pl[0].plot(sF.freqTime[freq],NP.rad2deg(clPhS[freq,:,0]%(2*NP.pi)),'.',markersize=3.,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[1].plot(sF.freqTime[freq],NP.rad2deg(clPhS[freq,:,3]%(2*NP.pi)),'.',markersize=3.,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[2].plot(sF.freqTime[freq],NP.rad2deg(clPhS[freq,:,6]%(2*NP.pi)),'.',markersize=3.,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
for pic in range(3):
    pl[pic].xaxis.set_major_locator(PL.MultipleLocator(120));
    pl[pic].xaxis.set_minor_locator(PL.MultipleLocator(30));
    pl[pic].xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
    pl[pic].set_xlim(1*3600,1.25*3600)
    pl[pic].set_ylim(0,400)
    if pic == 2:
        pl[pic].set_xlabel('UT')
    pl[pic].set_ylabel('degree')
    pl[pic].legend()
    pl[pic].grid()

fig, pl = PL.subplots(nrows=3,ncols=1,figsize=(8,8))
fig.tight_layout()
pl[0].set_title('N1-N2')
pl[0].text(sF.freqTime[0,1],360,'SRH phases %s'%(sF.dateObs.split('T')[0]),fontsize=14)
pl[1].set_title('N2-N3')
pl[2].set_title('N1-N3')
for freq in range(sF.freqListLength):
    pl[0].plot(sF.freqTime[freq],NP.rad2deg(NP.angle(sF.visLcp[freq,:,3007])%(2*NP.pi)),'.',markersize=3.,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[1].plot(sF.freqTime[freq],NP.rad2deg(NP.angle(sF.visLcp[freq,:,3007+1])%(2*NP.pi)),'.',markersize=3.,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[2].plot(sF.freqTime[freq],NP.rad2deg(NP.angle(sF.visLcp[freq,:,3007+30])%(2*NP.pi)),'.',markersize=3.,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
for pic in range(3):
    pl[pic].xaxis.set_major_locator(PL.MultipleLocator(120));
    pl[pic].xaxis.set_minor_locator(PL.MultipleLocator(30));
    pl[pic].xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
    pl[pic].set_xlim(1*3600,1.25*3600)
    pl[pic].set_ylim(0,400)
    if pic == 2:
        pl[pic].set_xlabel('UT')
    pl[pic].set_ylabel('degree')
    pl[pic].legend()
    pl[pic].grid()

fig, pl = PL.subplots(nrows=3,ncols=1,figsize=(8,8))
fig.tight_layout()
pl[0].set_title('N1-N2')
pl[0].text(sF.freqTime[0,1],0.09,'SRH amplitudes %s'%(sF.dateObs.split('T')[0]),fontsize=14)
pl[1].set_title('N2-N3')
pl[2].set_title('N1-N3')
for freq in range(sF.freqListLength):
    pl[0].plot(sF.freqTime[freq],NP.abs(sF.visLcp[freq,:,3007])%(2*NP.pi),'.',markersize=3.,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[1].plot(sF.freqTime[freq],NP.abs(sF.visLcp[freq,:,3007+1])%(2*NP.pi),'.',markersize=3.,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[2].plot(sF.freqTime[freq],NP.abs(sF.visLcp[freq,:,3007+30])%(2*NP.pi),'.',markersize=3.,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
for pic in range(3):
    pl[pic].xaxis.set_major_locator(PL.MultipleLocator(120));
    pl[pic].xaxis.set_minor_locator(PL.MultipleLocator(30));
    pl[pic].xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
    pl[pic].set_xlim(1*3600,1.25*3600)
    pl[pic].set_ylim(0,0.1)
    if pic == 2:
        pl[pic].set_xlabel('UT')
    pl[pic].set_ylabel('correlation')
    pl[pic].legend()
    pl[pic].grid()
