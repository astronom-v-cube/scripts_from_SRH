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
#dateTime = '20201031T004345'
#dateTime = '20201031T030324'
dateTime = '20201031T041418'

antNumber = 30
ssrtAnt = NP.linspace(0,antNumber-1,antNumber,dtype='int')*2+197
srhAnt = NP.linspace(0,antNumber-1,antNumber,dtype='int') + 1
amps_lcp = []
amps_rcp = []

for i in range(antNumber):
    amps_lcp.append(fits.open(('srh_amp_%d_' + dateTime + '.fits')%(ssrtAnt[i]))[1].data['amp_lcp']*ampScale)
    amps_rcp.append(fits.open(('srh_amp_%d_' + dateTime + '.fits')%(ssrtAnt[i]))[1].data['amp_rcp']*ampScale)
    
amp_211 = fits.open('srh_amp_211_' + dateTime + '.fits')

freqList = amp_211[1].data['frequency']
time = amp_211[1].data['time']
global _timeStart
_timeStart = time[0]
sampleNumber = time.shape[0]

grayColors = PL.get_cmap('gray')

PL.figure(figsize=(10,10))
pl0 = PL.subplot(111)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
for i in range(antNumber):
    if NP.array(amps_lcp[i]).mean() < 0.1:
        pl0.plot(time,amps_lcp[i], label='ant %d (N%02d)'% (ssrtAnt[i], srhAnt[i]), color = grayColors(int(i/(antNumber-1)*150)))
    else:
        if i == 19:
            pl0.plot(time, (amps_lcp[i] + amps_rcp[i])*.5, linewidth = 0.5, color = 'red')
        else:
            pl0.plot(time, (amps_lcp[i] + amps_rcp[i])*.5, linewidth = 0.5, color = grayColors(int(i/(antNumber-1)*150)))

pl0.set_title(dateTime + ', %d MHz' % int(freqList[0]))
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
meanI = (meanLcp + meanRcp)*.5
meanV = (meanLcp - meanRcp)*.5
meanI -= meanI.min()
sfuScale = 70 / meanI.max()
meanI *= sfuScale 
meanV *= sfuScale

PL.figure(figsize=(10,10))
pl0 = PL.subplot(111)
#pl0.set_ylim(0.42,0.46)
#pl0.set_ylim(-0.1,0.5)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl0.plot(time, meanI)
pl0.plot(time, meanV)
pl0.set_xlabel('UTC')
pl0.set_title(dateTime + ', %d MHz' % int(freqList[0]))
pl0.grid()
