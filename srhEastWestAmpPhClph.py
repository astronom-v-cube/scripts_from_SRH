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

fitNames = findFits('SRH36_temp_20210504_2','*.fit')
fitNames.sort();
fitNames = fitNames[1:];

visScale = 1/(2e6*49)
ampScale = visScale / 128
firstBaseline = 3472 + 30
secondBaselineOffset = 96

for fName in fitNames:
    if (fName == fitNames[0]):
        sF = SrhFitsFile(fName, 257)
    else:
        sF.append(fName)

T0 = .5*3600
T1 = 2.5*3600
titleX = T0 + 10
titleY = 300
markSize = 2.
clPhEW = NP.zeros((sF.freqListLength,sF.freqTime.shape[1],95))
for freq in range(sF.freqListLength):
    firstBaseline = 3472
    secondBaselineOffset = 96
    for i in range(95):
        phaseA1A2 = NP.angle(sF.visLcp[freq,:,firstBaseline + i])
        phaseA2A3 = NP.angle(sF.visLcp[freq,:,firstBaseline + 1 + i])
        phaseA1A3 = NP.angle(sF.visLcp[freq,:,firstBaseline + secondBaselineOffset + i])
        clPhEW[freq,:,i] = phaseA1A2 + phaseA2A3 - phaseA1A3
firstBaseline = 3472 + 30
#-------------------------------------------------------------------------------------------------------------------------------------
fig, pl = PL.subplots(nrows=3,ncols=1,figsize=(8,8))
fig.tight_layout()
pl[0].set_title('triplet W2+W1-C0')
pl[0].text(titleX,titleY,'SRH closure phases %s'%(sF.dateObs.split('T')[0]),fontsize=14)
pl[1].set_title('triplet W1+C0-E1')
pl[2].set_title('triplet C0+E1-E2')
for freq in range(sF.freqListLength):
    pl[0].plot(sF.freqTime[freq],NP.rad2deg(clPhEW[freq,:,30]%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[1].plot(sF.freqTime[freq],NP.rad2deg(clPhEW[freq,:,31]%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[2].plot(sF.freqTime[freq],NP.rad2deg(clPhEW[freq,:,32]%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
for pic in range(3):
    pl[pic].xaxis.set_major_locator(PL.MultipleLocator(600));
    pl[pic].xaxis.set_minor_locator(PL.MultipleLocator(60));
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
pl[0].set_title('triplet E1+E2-E3')
pl[0].text(titleX,titleY,'SRH closure phases %s'%(sF.dateObs.split('T')[0]),fontsize=14)
pl[1].set_title('triplet E2+E3-E4')
pl[2].set_title('triplet E3+E4-E5')
for freq in range(sF.freqListLength):
    pl[0].plot(sF.freqTime[freq],NP.rad2deg(clPhEW[freq,:,33]%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[1].plot(sF.freqTime[freq],NP.rad2deg(clPhEW[freq,:,34]%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[2].plot(sF.freqTime[freq],NP.rad2deg(clPhEW[freq,:,35]%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
for pic in range(3):
    pl[pic].xaxis.set_major_locator(PL.MultipleLocator(600));
    pl[pic].xaxis.set_minor_locator(PL.MultipleLocator(60));
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
pl[0].set_title('pair W2-W1')
pl[0].text(titleX,titleY,'SRH visibility phases %s'%(sF.dateObs.split('T')[0]),fontsize=14)
pl[1].set_title('pair W1-C0')
pl[2].set_title('pair C0-E1')
for freq in range(sF.freqListLength):
    pl[0].plot(sF.freqTime[freq],NP.rad2deg(NP.angle(sF.visLcp[freq,:,firstBaseline])%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[1].plot(sF.freqTime[freq],NP.rad2deg(NP.angle(sF.visLcp[freq,:,firstBaseline+1])%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[2].plot(sF.freqTime[freq],NP.rad2deg(NP.angle(sF.visLcp[freq,:,firstBaseline+2])%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
for pic in range(3):
    pl[pic].xaxis.set_major_locator(PL.MultipleLocator(600));
    pl[pic].xaxis.set_minor_locator(PL.MultipleLocator(60));
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
pl[0].set_title('pair E1-E2')
pl[0].text(titleX,titleY,'SRH visibility phases %s'%(sF.dateObs.split('T')[0]),fontsize=14)
pl[1].set_title('pair E2-E3')
pl[2].set_title('pair E3-E4')
for freq in range(sF.freqListLength):
    pl[0].plot(sF.freqTime[freq],NP.rad2deg(NP.angle(sF.visLcp[freq,:,firstBaseline+3])%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[1].plot(sF.freqTime[freq],NP.rad2deg(NP.angle(sF.visLcp[freq,:,firstBaseline+4])%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[2].plot(sF.freqTime[freq],NP.rad2deg(NP.angle(sF.visLcp[freq,:,firstBaseline+5])%(2*NP.pi)),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
for pic in range(3):
    pl[pic].xaxis.set_major_locator(PL.MultipleLocator(600));
    pl[pic].xaxis.set_minor_locator(PL.MultipleLocator(60));
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
pl[0].set_title('pair W2-W1')
pl[0].text(titleX,0.27,'SRH visibility amplitudes %s'%(sF.dateObs.split('T')[0]),fontsize=14)
pl[1].set_title('pair W1-C0')
pl[2].set_title('pair C0-E1')
for freq in range(sF.freqListLength):
    pl[0].plot(sF.freqTime[freq],NP.abs(sF.visLcp[freq,:,firstBaseline])%(2*NP.pi),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[1].plot(sF.freqTime[freq],NP.abs(sF.visLcp[freq,:,firstBaseline+1])%(2*NP.pi),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[2].plot(sF.freqTime[freq],NP.abs(sF.visLcp[freq,:,firstBaseline+2])%(2*NP.pi),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
for pic in range(3):
    pl[pic].xaxis.set_major_locator(PL.MultipleLocator(600));
    pl[pic].xaxis.set_minor_locator(PL.MultipleLocator(60));
    pl[pic].xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
    pl[pic].set_xlim(T0,T1)
    pl[pic].set_ylim(0,0.3)
    if pic == 2:
        pl[pic].set_xlabel('UT')
    pl[pic].set_ylabel('correlation')
    pl[pic].legend(markerscale=3)
    pl[pic].grid()

fig, pl = PL.subplots(nrows=3,ncols=1,figsize=(8,8))
fig.tight_layout()
pl[0].set_title('pair E1-E2')
pl[0].text(titleX,0.27,'SRH visibility amplitudes %s'%(sF.dateObs.split('T')[0]),fontsize=14)
pl[1].set_title('pair E2-E3')
pl[2].set_title('pair E3-E4')
for freq in range(sF.freqListLength):
    pl[0].plot(sF.freqTime[freq],NP.abs(sF.visLcp[freq,:,firstBaseline+3])%(2*NP.pi),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[1].plot(sF.freqTime[freq],NP.abs(sF.visLcp[freq,:,firstBaseline+4])%(2*NP.pi),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
    pl[2].plot(sF.freqTime[freq],NP.abs(sF.visLcp[freq,:,firstBaseline+5])%(2*NP.pi),'.',markersize=markSize,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
for pic in range(3):
    pl[pic].xaxis.set_major_locator(PL.MultipleLocator(600));
    pl[pic].xaxis.set_minor_locator(PL.MultipleLocator(60));
    pl[pic].xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
    pl[pic].set_xlim(T0,T1)
    pl[pic].set_ylim(0,0.3)
    if pic == 2:
        pl[pic].set_xlabel('UT')
    pl[pic].set_ylabel('correlation')
    pl[pic].legend(markerscale=3)
    pl[pic].grid()
#-------------------------------------------------------------------------------------------------------------------------------------

#fig, pl = PL.subplots(nrows=1,ncols=1,figsize=(8,8))
#fig.tight_layout()
#pl.text(0,titleY,'EastWest mean closure phases %s'%(sF.dateObs.split('T')[0]),fontsize=14)
#for freq in range(sF.freqListLength):
#    pl.plot(NP.rad2deg((clPhEW[freq]%(2*NP.pi)).mean(0)),'.',markersize=9.,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
#pl.set_ylim(0,360)
#pl.set_xlabel('baseline triplet')
#pl.set_ylabel('degree')
#pl.legend()
#pl.grid()
#
#fig, pl = PL.subplots(nrows=1,ncols=1,figsize=(8,8))
#fig.tight_layout()
#pl.text(0,0.09,'EastWest mean visibility amplitudes %s'%(sF.dateObs.split('T')[0]),fontsize=14)
#for freq in range(sF.freqListLength):
#    pl.plot(NP.abs(sF.visLcp[freq,:,firstBaseline:firstBaseline+96]).mean(0),'.',markersize=9.,label='%0.1f GHz'%(sF.freqList[freq]*1e-6))
#pl.set_ylim(0,0.3)
#pl.set_xlabel('pair')
#pl.set_ylabel('correlation')
#pl.legend()
#pl.grid()
