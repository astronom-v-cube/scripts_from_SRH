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
dateTime = '20201017T023545'

amp_211 = fits.open('srh_amp_211_' + dateTime + '.fits')
amp_213 = fits.open('srh_amp_213_' + dateTime + '.fits')
vis_211_213 = fits.open('srh_vis_211_213_' + dateTime + '.fits')

time = amp_211[1].data['time']
global _timeStart
_timeStart = time[0]
sampleNumber = time.shape[0]

amp_211_lcp = amp_211[1].data['amp_lcp']*ampScale
amp_213_lcp = amp_213[1].data['amp_lcp']*ampScale

amp_211_smoothed = NP.convolve(amp_211_lcp, lf_filter, mode='same') / LF_size
amp_213_smoothed = NP.convolve(amp_213_lcp, lf_filter, mode='same') / LF_size
amp_211_noise = amp_211_lcp[LF_size//2:-LF_size//2] - amp_211_smoothed[LF_size//2:-LF_size//2]
amp_213_noise = amp_213_lcp[LF_size//2:-LF_size//2] - amp_213_smoothed[LF_size//2:-LF_size//2]

vis_211_213_lcp = vis_211_213[1].data['vis_lcp'][0:sampleNumber]*visScale
vis_211_213_lcp_vf_real = NP.sin(NP.pi/2*vis_211_213_lcp.real)
vis_211_213_lcp_vf_imag = NP.sin(NP.pi/2*vis_211_213_lcp.imag)
#vis_amp_211_213_lcp = vis_211_213_lcp_vf_real
vis_amp_211_213_lcp = vis_211_213_lcp_vf_imag
#vis_amp_211_213_lcp = NP.sqrt(vis_211_213_lcp_vf_real**2 + vis_211_213_lcp_vf_imag**2)
vis_211_213_smoothed = NP.convolve(vis_amp_211_213_lcp, lf_filter, mode='same') / LF_size
vis_211_213_noise = vis_amp_211_213_lcp[LF_size//2:-LF_size//2] - vis_211_213_smoothed[LF_size//2:-LF_size//2]

time = time[LF_size//2:-LF_size//2]
amp_211_smoothed = amp_211_smoothed[LF_size//2:-LF_size//2]
amp_213_smoothed = amp_213_smoothed[LF_size//2:-LF_size//2]
vis_211_213_smoothed = vis_211_213_smoothed[LF_size//2:-LF_size//2]

vis_211_ind = [1000,1200]
vis_213_ind = [3000,3200]
vis_211_213_ind = [5900,6100]

SNR_211 = NP.abs(vis_211_213_smoothed[vis_211_ind[0]:vis_211_ind[1]].mean()/vis_211_213_noise[vis_211_ind[0]:vis_211_ind[1]].std())
SNR_213 = NP.abs(vis_211_213_smoothed[vis_213_ind[0]:vis_213_ind[1]].mean()/vis_211_213_noise[vis_213_ind[0]:vis_213_ind[1]].std())
SNR_211_213 = NP.abs(vis_211_213_smoothed[vis_211_213_ind[0]:vis_211_213_ind[1]].mean()/vis_211_213_noise[vis_211_213_ind[0]:vis_211_213_ind[1]].std())

PL.figure(figsize=(10,10))
pl0 = PL.subplot(111)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl0.plot(time,amp_211_smoothed,label='antenna 211')
pl0.plot(time,amp_213_smoothed,label='antenna 213')
pl0.plot(time,vis_211_213_smoothed,label='vis amp 211-213')
pl0.plot(time,vis_211_213_noise,label='vis noise 211-213')
pl0.plot(time[vis_211_ind[0]:vis_211_ind[1]],vis_211_213_noise[vis_211_ind[0]:vis_211_ind[1]],'.',label='SNR 211 %.1f' % SNR_211)
pl0.plot(time[vis_213_ind[0]:vis_213_ind[1]],vis_211_213_noise[vis_213_ind[0]:vis_213_ind[1]],'.',label='SNR 213 %.1f' % SNR_213)
pl0.plot(time[vis_211_213_ind[0]:vis_211_213_ind[1]],vis_211_213_noise[vis_211_213_ind[0]:vis_211_213_ind[1]],'.',label='SNR 211-213 %.1f' % SNR_211_213)

pl0.set_title(dateTime + ', 3000 MHz SNR %.3f,%.3f,%.3f' % (SNR_211, SNR_213, SNR_211_213))
pl0.set_xlabel('UTC')
pl0.grid()
pl0.legend()

singleNoise = vis_211_213_noise[200:200 + 2000]
pairNoise = vis_211_213_noise[5300:5300 + 2000]
#pairNoise = vis_211_213_lcp_vf_imag[5300:5300 + 2000]
freqs = (NP.linspace(0,1999,2000)-1000) * 1/(2000*.1)

PL.figure(figsize=(10,10))
pl0 = PL.subplot(111)
pl0.plot(freqs,NP.roll(NP.abs(NP.fft.fft(singleNoise)),1000))
pl0.plot(freqs,NP.roll(NP.abs(NP.fft.fft(pairNoise)),1000))
pl0.grid()
pl0.set_xlabel('Hz')
pl0.set_ylabel('amplitude')
pl0.set_title('Spectrum of correlation fluctuations')

