#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 08:15:50 2020

@author: svlesovoi
"""
import numpy as NP;
import pylab as PL;
from astropy.io import fits
from BadaryRAO import BadaryRAO

dt_major = 3600.;
dt_minor = 900.;

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

RAO = BadaryRAO('2020-04-23')

redVis = fits.open('srh_redVis_20200423.fits')
redFreq = redVis[1].data['frequencies']
redTime = redVis[2].data['time']
redA = redVis[3].data['antennaA']
redB = redVis[3].data['antennaB']
redRcp = redVis[4].data['VIS_RCP']
redLcp = redVis[4].data['VIS_LCP']
redRcpArr = redRcp.reshape(redFreq.shape[0],redTime.shape[1],redA.shape[0])
redLcpArr = redLcp.reshape(redFreq.shape[0],redTime.shape[1],redA.shape[0])

vis = 64
freq = 30

fig = PL.figure(1,figsize = (20,16));
fig.suptitle('Antenna pair %d, %d, frequency %.3f MHz'%(redA[vis], redB[vis], redFreq[freq]))
spl = fig.add_subplot(1,1,1)
spl.set_xlabel('UTC')
#spl.axes.xaxis.set_major_locator(PL.MultipleLocator(dt_major))
#spl.axes.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
#spl.axes.xaxis.set_minor_locator(PL.MultipleLocator(dt_minor))

#spl.plot(redTime[0,:], NP.abs(redRcpArr[freq,:,vis]))
#spl.plot(redTime[0,:], NP.abs(redLcpArr[freq,:,vis]))
#spl.plot(redTime[0,:], NP.angle(redRcpArr[freq,:,vis]), '.')
#spl.plot(redTime[0,:], NP.angle(redLcpArr[freq,:,vis]), '.')
#spl.plot(redTime[freq], redRcpArr[freq,:,vis].real)
#spl.plot(redTime[freq], redRcpArr[freq,:,vis].imag)
#spl.plot(redTime[freq], redLcpArr[freq,:,vis].real)
#spl.plot(redTime[freq], redLcpArr[freq,:,vis].imag)
#for f in range(redFreq.shape[0]):
for f in NP.linspace(25,30,2,dtype='int'):
    hA = NP.deg2rad(15*(redTime[f,:] - RAO.culmination)/3600)
    cosP = NP.sin(hA) * NP.cos(RAO.declination)
    theta = NP.rad2deg(3e8/(redFreq[f]*1e6)/(9.6*NP.sqrt(1 - cosP*cosP)))*60*NP.sign(hA)
    spl.plot(theta, NP.abs(redRcpArr[f,:,vis]) + NP.abs(redLcpArr[f,:,vis]),linewidth=0.3)
spl.grid()
