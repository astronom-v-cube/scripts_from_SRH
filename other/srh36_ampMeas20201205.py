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
#dateTime = '20201205T020236'
#dateTime = '20201205T030947'
dateTime = '20201205T064815'


antNumber = 29
ant = NP.linspace(0,antNumber-1,antNumber,dtype='int')*2+199
amps_lcp = []
amps_rcp = []

amp_0 = fits.open(('srh_amp_%d_' + dateTime + '.fits')%(ant[0]))

time = amp_0[1].data['time']
frequency = amp_0[1].data['frequency']
sampleNumber = time.shape[1]
amps_lcp = amp_0[1].data['amp_lcp']*ampScale
amps_rcp = amp_0[1].data['amp_rcp']*ampScale

for i in NP.linspace(1,28,28,dtype='int'):
    amp_i = fits.open(('srh_amp_%d_' + dateTime + '.fits')%(ant[i]))
    amps_lcp = NP.concatenate((amps_lcp, amp_i[1].data['amp_lcp']*ampScale),axis=1)
    amps_rcp = NP.concatenate((amps_rcp, amp_i[1].data['amp_rcp']*ampScale),axis=1)
    
amps_lcp = amps_lcp.reshape(frequency.shape[0], antNumber, sampleNumber)
amps_rcp = amps_rcp.reshape(frequency.shape[0], antNumber, sampleNumber)

PL.figure(figsize=(10,10))
pl0 = PL.subplot(111)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
for i in range(antNumber):
#    for j in range(frequency.shape[0]):
#        pl0.plot(time[j][1:], amps_lcp[j,i][1:] + amps_rcp[j,i][1:], color='red')
#        pl0.plot(time[j][1:], amps_lcp[j,i][1:] - amps_rcp[j,i][1:], color='blue')
     pl0.plot(time[0][1:], amps_lcp[0,i][1:] + amps_rcp[0,i][1:], color='red')

#pl0.set_title(dateTime + ', 3030 MHz')
pl0.set_xlabel('UTC')
pl0.grid()
#pl0.legend()

meanLcp = amps_lcp.mean(axis=1)
meanRcp = amps_rcp.mean(axis=1)

meanLInterp = NP.zeros((meanLcp.shape[0],meanLcp.shape[1]*2))
meanRInterp = NP.zeros((meanLcp.shape[0],meanLcp.shape[1]*2))
meanI = NP.zeros((meanLcp.shape[0],meanLcp.shape[1]*2))
meanV = NP.zeros((meanLcp.shape[0],meanLcp.shape[1]*2))
timeInterp = NP.linspace(time[0,0],time[0,-1],time.shape[1]*2)

flux = NP.linspace(70,120,frequency.shape[0])

for i in range(meanLcp.shape[0]):
    meanLInterp[i] = NP.interp(timeInterp,time[i],meanLcp[i] - 0.2)
    meanRInterp[i] = NP.roll(NP.interp(timeInterp,time[i],meanRcp[i] - 0.2),-1,axis=0)
    meanI[i] = (meanLInterp[i] + meanRInterp[i])*.5
    sfuScale = flux[i] / meanI[i][0:1000].max()
    meanI[i] *= sfuScale 
    meanV[i] = (meanLInterp[i] - meanRInterp[i])*.5 * sfuScale
    meanV[i] -= meanV[i].mean()


sens =   meanV[0].std()
PL.figure(figsize=(10,10))
pl0 = PL.subplot(111)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
for i in range(frequency.shape[0]):
    pl0.plot(timeInterp[:-2], meanI[i][:-2], label = 'I%d'%(frequency[i]/1000))
    pl0.plot(timeInterp[:-2], meanV[i][:-2], label = 'V%d'%(frequency[i]/1000))
#    pl0.plot(timeInterp[:-2], meanLInterp[i][:-2], label = 'L%d'%(i))
#    pl0.plot(timeInterp[:-2], meanRInterp[i][:-2], label = 'R%d'%(i))
pl0.set_ylabel('s.f.u.')
pl0.set_xlabel('UTC')
pl0.set_title(dateTime + ', %d MHz,'%frequency[0] + ' sensitivity = %.4f s.f.u.'%sens)
pl0.grid()
pl0.legend()

