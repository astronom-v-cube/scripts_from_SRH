#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 02:08:41 2021

@author: svlesovoi
"""

import json
import numpy as NP
import pylab as PL
from astropy.io import fits
from scipy import signal

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

win = signal.windows.gaussian(2199,2155)
rainbowColors = PL.get_cmap('rainbow')

#goesFile = open('xrays-6-hour.json')
goesFile = open('xrays-16-1-day.json')
#goesFile = open('xrays-1-day.json')
#goesFile = open('xrays-20210822.json')
goesData = json.load(goesFile)
goesFile.close()

x0822 = []
for x in goesData:
    if (x['time_tag'].split('T')[0] == '2021-09-12'): 
#    if (x['time_tag'].split('T')[0] == '2021-08-22'): 
#    if (x['time_tag'].split('T')[0] == '2021-08-13'): 
#    if (x['time_tag'].split('T')[0] == '2021-08-14'): 
#    if (x['time_tag'].split('T')[0] == '2021-08-15'): 
#    if (x['time_tag'].split('T')[0] == '2021-08-16'): 
#    if (x['time_tag'].split('T')[0] == '2021-08-17'): 
        x0822.append(x)
goesData = x0822
        
N = len(goesData)
xrays_time = NP.zeros(N//2)
xrays_4_8 = NP.zeros(N//2)
xrays_005_04 = NP.zeros(N//2)

for i in range(N):
    if i%2:
        xrays_4_8[i//2] = goesData[i]['flux']
        hhmm = goesData[i]['time_tag'].split('T')[1].split(':')[0:2]
        xrays_time[i//2] = 3600*int(hhmm[0]) + 60*int(hhmm[1])
#        if xrays_time[i//2] > 10*3600:
#            xrays_time[i//2] -= 24*3600
    else:
        xrays_005_04[i//2] = goesData[i]['flux']

#cF = fits.open('srh_cp_20210412.fits')

cF = fits.open('srh_cp_20210912.fits')
#cF = fits.open('srh_cp_20210822.fits')
#cF = fits.open('srh_cp_20210813.fits')
#cF = fits.open('srh_cp_20210814.fits')
#cF = fits.open('srh_cp_20210815.fits')
#cF = fits.open('srh_cp_20210816.fits')
#cF = fits.open('srh_cp_20210817.fits')

srhFreqList = cF[1].data['frequencies']
srhTime = cF[2].data['time']
srhCorrI = cF[2].data['I']
srhCorrV = cF[2].data['V']
srhFluxI = cF[2].data['flux_I']
srhMeanFluxI = srhFluxI.mean(axis=0)
srhMeanFluxISmoothed = signal.convolve(srhMeanFluxI,win,mode='same')/win.sum()

t0 = 0

fig = PL.figure(figsize=(12,7))
sub = fig.add_subplot(1,1,1);
sub.set_ylabel('flux');
sub.set_xlabel('UT');
sub.xaxis.set_major_locator(PL.MultipleLocator(3600));
sub.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
sub.xaxis.set_minor_locator(PL.MultipleLocator(600));
sub.set_xlim(0,10.0*3600)
sub.set_ylim(0,50e-7)

sub.plot(xrays_time[t0:],xrays_4_8[t0:],label='GOES X-Ray 0.1-0.8 nm',color='red',markersize=0.2)
sub.plot(xrays_time[t0:],xrays_005_04[t0:],label='GOES X-Ray 0.05-0.4 nm',color='blue',markersize=0.2)
#for freq in range(srhFreqList.shape[0]):
#    sub.plot(srhTime[freq],srhCorrI[freq]*1e-4,'.',markersize=0.2,color=rainbowColors(100+(srhFreqList.shape[0] - freq)*20),label='SRH %d MHz'%(srhFreqList[freq]*1e-3))
sub.plot([0,10*3600],[1e-7,1e-7], label='X-ray flare class A')

tsub = sub.twinx()
tsub.set_ylabel('correlation');
tsub.set_ylim(0,0.015)
tsub.plot(srhTime[0],srhCorrI[0:3].mean(axis=0),color='orange',label='SRH 2.8-3.4 GHz',linewidth=0.8)
tsub.plot(srhTime[3],srhCorrI[3:6].mean(axis=0),color='green',label='SRH 3.9-5.6 GHz',linewidth=0.8)

sub.grid()
sub.legend(loc='upper left')
tsub.legend()
sub.set_title('SRH and GOES lightcurves %s'%(cF[0].header['DATE-OBS']))
