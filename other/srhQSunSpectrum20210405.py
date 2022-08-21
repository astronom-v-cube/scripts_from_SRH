#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 09:28:23 2021

@author: sergeyvlesovoi
"""

from astropy.io import fits
import pylab as PL
import numpy as NP
from scipy import signal

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

#cF = fits.open('/home/sergeyvlesovoi/SRH36/corrPlot/srh_cp_20210404.fits')
cF = fits.open('srh_cp_20210411.fits')
#fF = fits.open('/home/sergeyvlesovoi/SRH36/fluxPlot/srh_flux_20210404.fits')

win = signal.windows.gaussian(31,15)

t0 = 21000
t1 = 25000
ampScale = 2e-6
ampOffset = -0.004
freqs = cF[1].data['frequencies']
times = cF[2].data['time'][:,t0:t1]
cpI = cF[2].data['I'][:,t0:t1]
flI = cF[2].data['flux_I'][:,t0:t1]*ampScale

for f in range(freqs.shape[0]):
    cpI[f] = signal.convolve(cpI[f],win,mode='same')/win.sum()
    cpI[f] -= cpI[f].min()
    
    flI[f] = signal.convolve(flI[f],win,mode='same')/win.sum()
    flI[f] -= flI[f].min()
    flI[f] += ampOffset

t1 = t1 - t0 - 1300
t0 = 500
times = times[:,t0:t1]
cpI = cpI[:,t0:t1]
flI = flI[:,t0:t1]
fig = PL.figure()
sub = fig.add_subplot(1,1,1);
sub.set_ylabel('a.u.');
sub.set_xlabel('UT');
sub.xaxis.set_major_locator(PL.MultipleLocator(3600));
sub.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
sub.xaxis.set_minor_locator(PL.MultipleLocator(900));

for f in range(freqs.shape[0]):
    sub.plot(times[f], cpI[f])
    sub.plot(times[f], flI[f])
sub.plot(times[f], flI.sum(0)/freqs.shape[0])
#sub.plot(times[f], flI.sum(0)/freqs.shape[0],'o')
sub.grid()


