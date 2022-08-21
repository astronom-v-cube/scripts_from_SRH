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
#dateTime = '20201012T010233'
#dateTime = '20201012T020911'
#dateTime = '20201012T025706'
dateTime = '20201012T034536'

amp_211 = fits.open('srh_amp_211_' + dateTime + '.fits')
amp_213 = fits.open('srh_amp_213_' + dateTime + '.fits')
#vis_211_213 = fits.open('srh_vis_211_213_' + dateTime + '.fits')

time = amp_211[1].data['time']
global _timeStart
_timeStart = time[0]
sampleNumber = time.shape[0]
#forestRange = [6350,6350+1200]
#forestRange = [2350,2350+1200]
forestRange = [1350,1350+200]
T_forest = 300

amp_211_lcp = amp_211[1].data['amp_lcp'] * ampScale
amp_213_lcp = amp_213[1].data['amp_lcp'] * ampScale
amp_211_lcp -= amp_211_lcp.min()
amp_211_lcp /= amp_211_lcp.max()
amp_213_lcp -= amp_213_lcp.min()
amp_213_lcp /= amp_213_lcp.max()

amp_211_smoothed = NP.convolve(amp_211_lcp, lf_filter, mode='same') / LF_size
amp_213_smoothed = NP.convolve(amp_213_lcp, lf_filter, mode='same') / LF_size
amp_211_noise = amp_211_lcp[LF_size//2:-LF_size//2] - amp_211_smoothed[LF_size//2:-LF_size//2]
amp_213_noise = amp_213_lcp[LF_size//2:-LF_size//2] - amp_213_smoothed[LF_size//2:-LF_size//2]

#vis_211_213_lcp = vis_211_213[1].data['vis_lcp'][0:sampleNumber] * visScale
#vis_211_213_lcp_vf_real = NP.sin(NP.pi/2*vis_211_213_lcp.real)
#vis_211_213_lcp_vf_imag = NP.sin(NP.pi/2*vis_211_213_lcp.imag)
#vis_amp_211_213_lcp = NP.sqrt(vis_211_213_lcp_vf_real**2 + vis_211_213_lcp_vf_imag**2)
#vis_211_213_smoothed = NP.convolve(vis_amp_211_213_lcp, lf_filter, mode='same') / LF_size
#vis_211_213_noise = vis_amp_211_213_lcp[LF_size//2:-LF_size//2] - vis_211_213_smoothed[LF_size//2:-LF_size//2]

time = time[LF_size//2:-LF_size//2]
amp_211_smoothed = amp_211_smoothed[LF_size//2:-LF_size//2]
amp_213_smoothed = amp_213_smoothed[LF_size//2:-LF_size//2]

T_211_scale = T_forest/amp_211_smoothed[forestRange[0]:forestRange[1]].mean()
T_213_scale = T_forest/amp_213_smoothed[forestRange[0]:forestRange[1]].mean()

PL.figure(figsize=(10,10))
pl0 = PL.subplot(111, )
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl0.plot(time,amp_211_smoothed*T_211_scale,label='antenna 211')
#pl0.plot(time[forestRange[0]:forestRange[1]],amp_211_smoothed[forestRange[0]:forestRange[1]],label='antenna 211')
pl0.plot(time,amp_213_smoothed*T_213_scale,label='antenna 213')
pl0.plot([time[0],time[-1]],[T_forest,T_forest], label='300 K')

pl0.set_ylim((0,3000.))
#pl0.set_title(dateTime + ', 3000 MHz')
#pl0.set_title(dateTime + ', 4500 MHz')
pl0.set_title(dateTime + ', 6000 MHz')
pl0.set_xlabel('UTC')
pl0.set_ylabel('antenna temperature')
pl0.grid()
pl0.legend()
#PL.savefig(dateTime + '_4500_211_213.png')
PL.savefig(dateTime + '_6000_211_213.png')
