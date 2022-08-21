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

LF_size = 13
lf_filter = NP.ones(LF_size)

visScale = 1/(2e6*49)
ampScale = visScale / 128
LF_size = 19
lf_filter = NP.ones(LF_size)
dateTime = '20201016T071742'

amp_215 = fits.open('srh_amp_215_' + dateTime + '.fits')
amp_217 = fits.open('srh_amp_217_' + dateTime + '.fits')
vis_215_217 = fits.open('srh_vis_215_217_' + dateTime + '.fits')

time = amp_215[1].data['time']
global _timeStart
_timeStart = time[0]
sampleNumber = time.shape[0]

amp_215_lcp = amp_215[1].data['amp_lcp']
amp_217_lcp = amp_217[1].data['amp_lcp']
#amp_215_lcp = amp_215_lcp - amp_215_lcp.min()
#amp_215_lcp = amp_215_lcp / amp_215_lcp.max()
#amp_217_lcp = amp_217_lcp - amp_217_lcp.min()
#amp_217_lcp = amp_217_lcp / amp_217_lcp.max()

amp_215_smoothed = NP.convolve(amp_215_lcp, lf_filter, mode='same') / LF_size
amp_217_smoothed = NP.convolve(amp_217_lcp, lf_filter, mode='same') / LF_size
amp_215_noise = amp_215_lcp[LF_size//2:-LF_size//2] - amp_215_smoothed[LF_size//2:-LF_size//2]
amp_217_noise = amp_217_lcp[LF_size//2:-LF_size//2] - amp_217_smoothed[LF_size//2:-LF_size//2]

vis_215_217_lcp = vis_215_217[1].data['vis_lcp'][0:sampleNumber]
vis_215_217_lcp_vf_real = NP.sin(NP.pi/2*vis_215_217_lcp.real)
vis_215_217_lcp_vf_imag = NP.sin(NP.pi/2*vis_215_217_lcp.imag)
vis_amp_215_217_lcp = NP.sqrt(vis_215_217_lcp_vf_real**2 + vis_215_217_lcp_vf_imag**2)
vis_215_217_smoothed = NP.convolve(vis_amp_215_217_lcp, lf_filter, mode='same') / LF_size
vis_215_217_noise = vis_amp_215_217_lcp[LF_size//2:-LF_size//2] - vis_215_217_smoothed[LF_size//2:-LF_size//2]

time = time[LF_size//2:-LF_size//2]
amp_215_smoothed = amp_215_smoothed[LF_size//2:-LF_size//2]
amp_217_smoothed = amp_217_smoothed[LF_size//2:-LF_size//2]
vis_215_217_smoothed = vis_215_217_smoothed[LF_size//2:-LF_size//2]

PL.figure(figsize=(10,10))
pl0 = PL.subplot(111, )
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl0.plot(time,amp_215_smoothed,label='antenna 215')
pl0.plot(time,amp_217_smoothed,label='antenna 217')
pl0.plot(time,vis_215_217_smoothed,label='vis 215-217')

#pl0.set_ylim((0,3000.))
pl0.set_title(dateTime + ', 3000 MHz')
pl0.set_xlabel('UTC')
pl0.grid()
pl0.legend()
