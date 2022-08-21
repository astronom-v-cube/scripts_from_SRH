#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 01:36:25 2020

@author: svlesovoi
"""
import numpy as NP
import pylab as PL
from astropy.io import fits

def hhmm_format(t, pos):
  t = t#0.1*t + _timeStart
  hh = int(t / 3600.)
  t -= hh*3600.
  mm = int(t / 60.)
  t -= mm*60
  ss = int(t)
  return '%02d:%02d:%02d' % (hh,mm,ss)

def antenna_format(ant, pos):
    return '%03d'%(ant*2 + 199)

LF_size = 13
lf_filter = NP.ones(LF_size)

visScale = 1/(2e6*49)
ampScale = visScale / 128
#dateTime = '20201202T024031'
dateTime = '20201202T034032'


antNumber = 29
ant = NP.linspace(0,antNumber-1,antNumber)*2+199
amps_lcp = []
amps_rcp = []

for i in range(antNumber):
    amps_lcp.append(fits.open(('srh_amp_%d_' + dateTime + '.fits')%(ant[i]))[1].data['amp_lcp']*ampScale)
    amps_rcp.append(fits.open(('srh_amp_%d_' + dateTime + '.fits')%(ant[i]))[1].data['amp_rcp']*ampScale)
    
amp_211 = fits.open('srh_amp_211_' + dateTime + '.fits')

time = amp_211[1].data['time']
frequency = amp_211[1].data['frequency']
global _timeStart
_timeStart = time[0]
sampleNumber = time.shape[0]


PL.figure(figsize=(10,10))
pl0 = PL.subplot(111)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
for i in range(antNumber):
    if NP.array(amps_lcp[i]).mean() < 0.01:
        pl0.plot(time,amps_lcp[i], label='ant %d'%ant[i])
    else:
        pl0.plot(time, amps_lcp[i], color='red')
        pl0.plot(time, amps_rcp[i], color='blue')

pl0.set_title(dateTime + ', 3030 MHz')
pl0.set_xlabel('UTC')
pl0.grid()
pl0.legend()

#PL.figure()
#pl1 = PL.subplot(111)
#pl1.yaxis.set_major_formatter(PL.FuncFormatter(antenna_format))
#pl1.imshow(amps_lcp, aspect=100)
#
meanLcp = NP.array(amps_lcp).mean(0)
meanRcp = NP.array(amps_rcp).mean(0)

meanL = meanLcp - 0.2# 0.2 sky level from prevous calibration
meanR = meanRcp - 0.2

timeInterp = NP.linspace(time[0],time[-1],time.shape[0]*2)
meanLInterp = NP.interp(timeInterp,time,meanL)
meanRInterp = NP.roll(NP.interp(timeInterp,time,meanR),-1)

meanI = (meanLInterp + meanRInterp)*.5
sfuScale = 70 / meanI[0:1000].max()
meanI *= sfuScale 
meanV = (meanLInterp - meanRInterp)*.5 * sfuScale
meanV -= meanV.mean()
sens =   meanV.std()
meanL *= sfuScale
meanR *= sfuScale

PL.figure(figsize=(10,10))
pl0 = PL.subplot(111)
#pl0.set_ylim(0.42,0.46)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl0.plot(timeInterp[:-2], meanI[:-2], label = 'I')
pl0.plot(timeInterp[:-2], meanV[:-2], label = 'V')
#pl0.plot(timeInterp[:-2], meanLInterp[:-2], label = 'L')
#pl0.plot(timeInterp[:-2], meanRInterp[:-2], label = 'R')
pl0.set_ylabel('s.f.u.')
pl0.set_xlabel('UTC')
pl0.set_title(dateTime + ', %d MHz,'%frequency[0] + ' sensitivity = %.4f s.f.u.'%sens)
pl0.grid()
pl0.legend()

