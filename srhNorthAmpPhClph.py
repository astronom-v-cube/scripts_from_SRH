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

T0 = 1.03*3600
T1 = 1.73*3600
titleX = T0 + 10
titleY = 300
markSize = 2.
clPhS = NP.zeros((sF.freqListLength,sF.freqTime.shape[1],29))
clPhEW = NP.zeros((sF.freqListLength,sF.freqTime.shape[1],95))
for freq in range(sF.freqListLength):
    firstBaseline = 3007
    secondBaselineOffset = 30
    for i in range(29):
        phaseA1A2 = NP.angle(sF.visLcp[freq,:,firstBaseline + i])
        phaseA2A3 = NP.angle(sF.visLcp[freq,:,firstBaseline + 1 + i])
        phaseA1A3 = NP.angle(sF.visLcp[freq,:,firstBaseline + secondBaselineOffset + i])
        clPhS[freq,:,i] = phaseA1A2 + phaseA2A3 - phaseA1A3
    
    firstBaseline = 3472
    secondBaselineOffset = 96
    for i in range(95):
        phaseA1A2 = NP.angle(sF.visLcp[freq,:,firstBaseline + i])
        phaseA2A3 = NP.angle(sF.visLcp[freq,:,firstBaseline + 1 + i])
        phaseA1A3 = NP.angle(sF.visLcp[freq,:,firstBaseline + secondBaselineOffset + i])
        clPhEW[freq,:,i] = phaseA1A2 + phaseA2A3 - phaseA1A3
#-------------------------------------------------------------------------------------------------------------------------------------
fig, pl = PL.subplots(nrows=3,ncols=1,figsize=(8,8))
fig.tight_layout()
pl[0].set_title('triplet N1+N2-N3')
pl[0].text(titleX,titleY,'SRH closure phases %s'%(sF.dateObs.split('T')[0]),fontsize=14)
pl[1].set_title('triplet N2+N3-N4')
pl[2].set_title('triplet N3+N4-N5')
for freq in range(sF.freqListLength):
    pl[0].plot(sF.freqTime[freq],NP.rad2deg(clPhS[freq,:,0]%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[1].plot(sF.freqTime[freq],NP.rad2deg(clPhS[freq,:,1]%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[2].plot(sF.freqTime[freq],NP.rad2deg(clPhS[freq,:,2]%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
for pic in range(3):
    pl[pic].xaxis.set_major_locator(PL.MultipleLocator(120));
    pl[pic].xaxis.set_minor_locator(PL.MultipleLocator(30));
    pl[pic].xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
    pl[pic].yaxis.set_major_locator(PL.MultipleLocator(50));
    pl[pic].set_xlim(T0,T1)
    pl[pic].set_ylim(0,360)
    if pic == 2:
        pl[pic].set_xlabel('UT')
    pl[pic].set_ylabel('degree')
    pl[pic].legend(markerscale=3)
    pl[pic].grid()

fig, pl = PL.subplots(nrows=3,ncols=1,figsize=(8,8))
fig.tight_layout()
pl[0].set_title('triplet N4+N5-N6')
pl[0].text(titleX,titleY,'SRH closure phases %s'%(sF.dateObs.split('T')[0]),fontsize=14)
pl[1].set_title('triplet N5+N6-N7')
pl[2].set_title('triplet N6+N7-N8')
for freq in range(sF.freqListLength):
    pl[0].plot(sF.freqTime[freq],NP.rad2deg(clPhS[freq,:,3]%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[1].plot(sF.freqTime[freq],NP.rad2deg(clPhS[freq,:,4]%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[2].plot(sF.freqTime[freq],NP.rad2deg(clPhS[freq,:,5]%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
for pic in range(3):
    pl[pic].xaxis.set_major_locator(PL.MultipleLocator(120));
    pl[pic].xaxis.set_minor_locator(PL.MultipleLocator(30));
    pl[pic].xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
    pl[pic].yaxis.set_major_locator(PL.MultipleLocator(50));
    pl[pic].set_xlim(T0,T1)
    pl[pic].set_ylim(0,360)
    if pic == 2:
        pl[pic].set_xlabel('UT')
    pl[pic].set_ylabel('degree')
    pl[pic].legend(markerscale=3)
    pl[pic].grid()
#-------------------------------------------------------------------------------------------------------------------------------------

fig, pl = PL.subplots(nrows=3,ncols=1,figsize=(8,8))
fig.tight_layout()
pl[0].set_title('pair N1-N2')
pl[0].text(titleX,titleY,'SRH visibility phases %s'%(sF.dateObs.split('T')[0]),fontsize=14)
pl[1].set_title('pair N2-N3')
pl[2].set_title('pair N3-N4')
for freq in range(sF.freqListLength):
    pl[0].plot(sF.freqTime[freq],NP.rad2deg(NP.angle(sF.visLcp[freq,:,3007])%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[1].plot(sF.freqTime[freq],NP.rad2deg(NP.angle(sF.visLcp[freq,:,3007+1])%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[2].plot(sF.freqTime[freq],NP.rad2deg(NP.angle(sF.visLcp[freq,:,3007+2])%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
for pic in range(3):
    pl[pic].xaxis.set_major_locator(PL.MultipleLocator(120));
    pl[pic].xaxis.set_minor_locator(PL.MultipleLocator(30));
    pl[pic].xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
    pl[pic].yaxis.set_major_locator(PL.MultipleLocator(50));
    pl[pic].set_xlim(T0,T1)
    pl[pic].set_ylim(0,360)
    if pic == 2:
        pl[pic].set_xlabel('UT')
    pl[pic].set_ylabel('degree')
    pl[pic].legend(markerscale=3)
    pl[pic].grid()

fig, pl = PL.subplots(nrows=3,ncols=1,figsize=(8,8))
fig.tight_layout()
pl[0].set_title('pair N4-N5')
pl[0].text(titleX,titleY,'SRH visibility phases %s'%(sF.dateObs.split('T')[0]),fontsize=14)
pl[1].set_title('pair N5-N6')
pl[2].set_title('pair N6-N7')
for freq in range(sF.freqListLength):
    pl[0].plot(sF.freqTime[freq],NP.rad2deg(NP.angle(sF.visLcp[freq,:,3007+3])%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[1].plot(sF.freqTime[freq],NP.rad2deg(NP.angle(sF.visLcp[freq,:,3007+4])%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[2].plot(sF.freqTime[freq],NP.rad2deg(NP.angle(sF.visLcp[freq,:,3007+5])%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
for pic in range(3):
    pl[pic].xaxis.set_major_locator(PL.MultipleLocator(120));
    pl[pic].xaxis.set_minor_locator(PL.MultipleLocator(30));
    pl[pic].xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
    pl[pic].yaxis.set_major_locator(PL.MultipleLocator(50));
    pl[pic].set_xlim(T0,T1)
    pl[pic].set_ylim(0,360)
    if pic == 2:
        pl[pic].set_xlabel('UT')
    pl[pic].set_ylabel('degree')
    pl[pic].legend(markerscale=3)
    pl[pic].grid()
#-------------------------------------------------------------------------------------------------------------------------------------

fig, pl = PL.subplots(nrows=3,ncols=1,figsize=(8,8))
fig.tight_layout()
pl[0].set_title('pair N1-N2')
pl[0].text(titleX,0.09,'SRH visibility amplitudes %s'%(sF.dateObs.split('T')[0]),fontsize=14)
pl[1].set_title('pair N2-N3')
pl[2].set_title('pair N3-N4')
for freq in range(sF.freqListLength):
    pl[0].plot(sF.freqTime[freq],NP.abs(sF.visLcp[freq,:,3007])%(2*NP.pi),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[1].plot(sF.freqTime[freq],NP.abs(sF.visLcp[freq,:,3007+1])%(2*NP.pi),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[2].plot(sF.freqTime[freq],NP.abs(sF.visLcp[freq,:,3007+2])%(2*NP.pi),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
for pic in range(3):
    pl[pic].xaxis.set_major_locator(PL.MultipleLocator(120));
    pl[pic].xaxis.set_minor_locator(PL.MultipleLocator(30));
    pl[pic].xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
    pl[pic].set_xlim(T0,T1)
    pl[pic].set_ylim(0,0.1)
    if pic == 2:
        pl[pic].set_xlabel('UT')
    pl[pic].set_ylabel('correlation')
    pl[pic].legend(markerscale=3)
    pl[pic].grid()

fig, pl = PL.subplots(nrows=3,ncols=1,figsize=(8,8))
fig.tight_layout()
pl[0].set_title('pair N4-N5')
pl[0].text(titleX,0.09,'SRH visibility amplitudes %s'%(sF.dateObs.split('T')[0]),fontsize=14)
pl[1].set_title('pair N5-N6')
pl[2].set_title('pair N7-N8')
for freq in range(sF.freqListLength):
    pl[0].plot(sF.freqTime[freq],NP.abs(sF.visLcp[freq,:,3007+3])%(2*NP.pi),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[1].plot(sF.freqTime[freq],NP.abs(sF.visLcp[freq,:,3007+4])%(2*NP.pi),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[2].plot(sF.freqTime[freq],NP.abs(sF.visLcp[freq,:,3007+5])%(2*NP.pi),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
for pic in range(3):
    pl[pic].xaxis.set_major_locator(PL.MultipleLocator(120));
    pl[pic].xaxis.set_minor_locator(PL.MultipleLocator(30));
    pl[pic].xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
    pl[pic].set_xlim(T0,T1)
    pl[pic].set_ylim(0,0.1)
    if pic == 2:
        pl[pic].set_xlabel('UT')
    pl[pic].set_ylabel('correlation')
    pl[pic].legend(markerscale=3)
    pl[pic].grid()
#-------------------------------------------------------------------------------------------------------------------------------------

fig, pl = PL.subplots(nrows=1,ncols=1,figsize=(8,8))
fig.tight_layout()
pl.text(0,titleY,'North mean closure phases %s'%(sF.dateObs.split('T')[0]),fontsize=14)
for freq in range(sF.freqListLength):
    pl.plot(NP.rad2deg((clPhS[freq]%(2*NP.pi)).mean(0)),'o',markersize=9.,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl.plot([10,11,12],NP.rad2deg((clPhS[freq,:,10:13]%(2*NP.pi)).mean(0)),'x',markersize=13.,color='black')
    pl.plot([14,15,16],NP.rad2deg((clPhS[freq,:,14:17]%(2*NP.pi)).mean(0)),'x',markersize=13.,color='black')
    pl.plot([26,27,28],NP.rad2deg((clPhS[freq,:,26:]%(2*NP.pi)).mean(0)),'x',markersize=13.,color='black')
pl.set_ylim(0,360)
pl.set_xlabel('baseline triplet')
pl.set_ylabel('degree')
pl.legend()
pl.grid()

fig, pl = PL.subplots(nrows=1,ncols=1,figsize=(8,8))
fig.tight_layout()
pl.text(0,0.09,'North mean visibility amplitudes %s'%(sF.dateObs.split('T')[0]),fontsize=14)
for freq in range(sF.freqListLength):
    pl.plot(NP.abs(sF.visLcp[freq,:,3007:3007+30]).mean(0),'o',markersize=9.,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl.plot([11,12],NP.abs(sF.visLcp[freq,:,3007+11:3007+13]).mean(0),'x',markersize=13.,color='black')
    pl.plot([15,16],NP.abs(sF.visLcp[freq,:,3007+15:3007+17]).mean(0),'x',markersize=13.,color='black')
    pl.plot([27,28],NP.abs(sF.visLcp[freq,:,3007+27:3007+29]).mean(0),'x',markersize=13.,color='black')
pl.set_ylim(0,0.1)
pl.set_xlabel('pair')
pl.set_ylabel('correlation')
pl.legend()
pl.grid()
