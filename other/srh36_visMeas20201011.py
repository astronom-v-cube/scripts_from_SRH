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
dateTime = '20201011T003609'

amp_211 = fits.open('srh_amp_211_' + dateTime + '.fits')
amp_213 = fits.open('srh_amp_213_' + dateTime + '.fits')
vis_211_213 = fits.open('srh_vis_211_213_' + dateTime + '.fits')

amp_81 = fits.open('srh_amp_81_' + dateTime + '.fits')
amp_83 = fits.open('srh_amp_83_' + dateTime + '.fits')
vis_81_83 = fits.open('srh_vis_81_83_' + dateTime + '.fits')

amp_46 = fits.open('srh_amp_46_' + dateTime + '.fits')
amp_48 = fits.open('srh_amp_48_' + dateTime + '.fits')
vis_46_48 = fits.open('srh_vis_46_48_' + dateTime + '.fits')

amp_101 = fits.open('srh_amp_101_' + dateTime + '.fits')
amp_104 = fits.open('srh_amp_104_' + dateTime + '.fits')
vis_101_104 = fits.open('srh_vis_101_104_' + dateTime + '.fits')

time = amp_211[1].data['time']
global _timeStart
_timeStart = time[0]
sampleNumber = time.shape[0]
noiseRange = 150

amp_211_lcp = amp_211[1].data['amp_lcp'] * ampScale
amp_213_lcp = amp_213[1].data['amp_lcp'] * ampScale
amp_211_smoothed = NP.convolve(amp_211_lcp, lf_filter, mode='same') / LF_size
amp_213_smoothed = NP.convolve(amp_213_lcp, lf_filter, mode='same') / LF_size
amp_211_noise = amp_211_lcp[LF_size//2:-LF_size//2] - amp_211_smoothed[LF_size//2:-LF_size//2]
amp_213_noise = amp_213_lcp[LF_size//2:-LF_size//2] - amp_213_smoothed[LF_size//2:-LF_size//2]

vis_211_213_lcp = vis_211_213[1].data['vis_lcp'][0:sampleNumber] * visScale
vis_211_213_lcp_vf_real = NP.sin(NP.pi/2*vis_211_213_lcp.real)
vis_211_213_lcp_vf_imag = NP.sin(NP.pi/2*vis_211_213_lcp.imag)
vis_amp_211_213_lcp = NP.sqrt(vis_211_213_lcp_vf_real**2 + vis_211_213_lcp_vf_imag**2)
vis_211_213_smoothed = NP.convolve(vis_amp_211_213_lcp, lf_filter, mode='same') / LF_size
vis_211_213_noise = vis_amp_211_213_lcp[LF_size//2:-LF_size//2] - vis_211_213_smoothed[LF_size//2:-LF_size//2]

amp_81_lcp = amp_81[1].data['amp_lcp'] * ampScale
amp_83_lcp = amp_83[1].data['amp_lcp'] * ampScale
amp_81_lcp -= amp_81_lcp.min()
amp_81_lcp /= amp_81_lcp.max()
amp_83_lcp -= amp_83_lcp.min()
amp_83_lcp /= amp_83_lcp.max()
amp_81_smoothed = NP.convolve(amp_81_lcp, lf_filter, mode='same') / LF_size
amp_83_smoothed = NP.convolve(amp_83_lcp, lf_filter, mode='same') / LF_size
amp_81_noise = amp_81_lcp[LF_size//2:-LF_size//2] - amp_81_smoothed[LF_size//2:-LF_size//2]
amp_83_noise = amp_83_lcp[LF_size//2:-LF_size//2] - amp_83_smoothed[LF_size//2:-LF_size//2]

vis_81_83_lcp = vis_81_83[1].data['vis_lcp'][0:sampleNumber] * visScale
vis_81_83_lcp_vf_real = NP.sin(NP.pi/2*vis_81_83_lcp.real)
vis_81_83_lcp_vf_imag = NP.sin(NP.pi/2*vis_81_83_lcp.imag)
vis_amp_81_83_lcp = NP.sqrt(vis_81_83_lcp_vf_real**2 + vis_81_83_lcp_vf_imag**2)
vis_amp_81_83_lcp -= vis_amp_81_83_lcp.min()
vis_amp_81_83_lcp /= vis_amp_81_83_lcp.max()
vis_81_83_smoothed = NP.convolve(vis_amp_81_83_lcp, lf_filter, mode='same') / LF_size
vis_81_83_noise = vis_amp_81_83_lcp[LF_size//2:-LF_size//2] - vis_81_83_smoothed[LF_size//2:-LF_size//2]

amp_46_lcp = amp_46[1].data['amp_lcp'] * ampScale
amp_48_lcp = amp_48[1].data['amp_lcp'] * ampScale
amp_46_lcp -= amp_46_lcp.min()
amp_46_lcp /= amp_46_lcp.max()
amp_48_lcp -= amp_48_lcp.min()
amp_48_lcp /= amp_48_lcp.max()
amp_46_smoothed = NP.convolve(amp_46_lcp, lf_filter, mode='same') / LF_size
amp_48_smoothed = NP.convolve(amp_48_lcp, lf_filter, mode='same') / LF_size
amp_46_noise = amp_46_lcp[LF_size//2:-LF_size//2] - amp_46_smoothed[LF_size//2:-LF_size//2]
amp_48_noise = amp_48_lcp[LF_size//2:-LF_size//2] - amp_48_smoothed[LF_size//2:-LF_size//2]

vis_46_48_lcp = vis_46_48[1].data['vis_lcp'][0:sampleNumber] * visScale
vis_46_48_lcp_vf_real = NP.sin(NP.pi/2*vis_46_48_lcp.real)
vis_46_48_lcp_vf_imag = NP.sin(NP.pi/2*vis_46_48_lcp.imag)
vis_amp_46_48_lcp = NP.sqrt(vis_46_48_lcp_vf_real**2 + vis_46_48_lcp_vf_imag**2)
vis_amp_46_48_lcp -= vis_amp_46_48_lcp.min()
vis_amp_46_48_lcp /= vis_amp_46_48_lcp.max()
vis_46_48_smoothed = NP.convolve(vis_amp_46_48_lcp, lf_filter, mode='same') / LF_size
vis_46_48_noise = vis_amp_46_48_lcp[LF_size//2:-LF_size//2] - vis_46_48_smoothed[LF_size//2:-LF_size//2]

amp_101_lcp = amp_101[1].data['amp_lcp'] * ampScale
amp_104_lcp = amp_104[1].data['amp_lcp'] * ampScale
amp_101_lcp -= amp_101_lcp.min()
amp_101_lcp /= amp_101_lcp.max()
amp_104_lcp -= amp_104_lcp.min()
amp_104_lcp /= amp_104_lcp.max()
amp_101_smoothed = NP.convolve(amp_101_lcp, lf_filter, mode='same') / LF_size
amp_104_smoothed = NP.convolve(amp_104_lcp, lf_filter, mode='same') / LF_size
amp_101_noise = amp_101_lcp[LF_size//2:-LF_size//2] - amp_101_smoothed[LF_size//2:-LF_size//2]
amp_104_noise = amp_104_lcp[LF_size//2:-LF_size//2] - amp_104_smoothed[LF_size//2:-LF_size//2]

vis_101_104_lcp = vis_101_104[1].data['vis_lcp'][0:sampleNumber] * visScale
vis_101_104_lcp_vf_real = NP.sin(NP.pi/2*vis_101_104_lcp.real)
vis_101_104_lcp_vf_imag = NP.sin(NP.pi/2*vis_101_104_lcp.imag)
vis_amp_101_104_lcp = NP.sqrt(vis_101_104_lcp_vf_real**2 + vis_101_104_lcp_vf_imag**2)
vis_amp_101_104_lcp -= vis_amp_101_104_lcp.min()
vis_amp_101_104_lcp /= vis_amp_101_104_lcp.max()
vis_101_104_smoothed = NP.convolve(vis_amp_101_104_lcp, lf_filter, mode='same') / LF_size
vis_101_104_noise = vis_amp_101_104_lcp[LF_size//2:-LF_size//2] - vis_101_104_smoothed[LF_size//2:-LF_size//2]

PL.figure(figsize=(10,10))
pl0 = PL.subplot(111, )
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl0.plot(time[LF_size//2:-LF_size//2],vis_211_213_smoothed[LF_size//2:-LF_size//2], label='visibility 211-213 amp')
pl0.plot(time[LF_size//2:-LF_size//2],vis_211_213_noise, label='visibility 211-213 noise')
pl0.plot(time[LF_size//2:-LF_size//2],amp_211_smoothed[LF_size//2:-LF_size//2],label='antenna 211')
pl0.plot(time[LF_size//2:-LF_size//2],amp_213_smoothed[LF_size//2:-LF_size//2],label='antenna 213')
pl0.plot(time[LF_size//2:-LF_size//2],amp_211_noise,label='noise 211')
pl0.plot(time[LF_size//2:-LF_size//2],amp_213_noise,label='noise 213')

pl0.plot(time[LF_size//2:-LF_size//2-noiseRange],amp_211_noise[:-noiseRange])
pl0.plot(time[LF_size//2:-LF_size//2-noiseRange],amp_213_noise[:-noiseRange])

pl0.set_ylim((-.5,1.))
pl0.set_title(dateTime + ', 4500 MHz, amp SNR %.3f, vis SNR %.3f'%(amp_211_smoothed[:-noiseRange].mean()/amp_211_noise[:-noiseRange].std(), vis_211_213_smoothed[:-noiseRange].mean()/vis_211_213_noise[:-noiseRange].std()))
pl0.set_xlabel('UTC')
pl0.set_ylabel('correlation')
pl0.grid()
pl0.legend()
PL.savefig(dateTime + '_4500_211_213.png')

PL.figure(figsize=(10,10))
pl1 = PL.subplot(111)
pl1.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl1.plot(time[LF_size//2:-LF_size//2],vis_81_83_smoothed[LF_size//2:-LF_size//2], label='visibility 81-83 amp')
pl1.plot(time[LF_size//2:-LF_size//2],vis_81_83_noise, label='visibility 81-83 noise')
pl1.plot(time[LF_size//2:-LF_size//2],amp_81_smoothed[LF_size//2:-LF_size//2],label='antenna 81')
pl1.plot(time[LF_size//2:-LF_size//2],amp_83_smoothed[LF_size//2:-LF_size//2],label='antenna 83')
pl1.plot(time[LF_size//2:-LF_size//2],amp_81_noise,label='noise 81')
pl1.plot(time[LF_size//2:-LF_size//2],amp_83_noise,label='noise 83')

pl1.plot(time[LF_size//2:-LF_size//2-noiseRange],amp_81_noise[:-noiseRange])
pl1.plot(time[LF_size//2:-LF_size//2-noiseRange],amp_83_noise[:-noiseRange])

pl1.set_ylim((-.5,1.))
pl1.set_title(dateTime + ', 4500 MHz, amp SNR %.3f, vis SNR %.3f'%(amp_81_smoothed[:-noiseRange].mean()/amp_81_noise[:-noiseRange].std(), vis_81_83_smoothed[:-noiseRange].mean()/vis_81_83_noise[:-noiseRange].std()))
pl1.set_xlabel('UTC')
pl1.set_ylabel('correlation')
pl1.grid()
pl1.legend()
PL.savefig(dateTime + '_4500_81_83.png')

PL.figure(figsize=(10,10))
pl1 = PL.subplot(111)
pl1.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl1.plot(time[LF_size//2:-LF_size//2],vis_46_48_smoothed[LF_size//2:-LF_size//2], label='visibility 46-48 amp')
pl1.plot(time[LF_size//2:-LF_size//2],vis_46_48_noise, label='visibility 46-48 noise')
pl1.plot(time[LF_size//2:-LF_size//2],amp_46_smoothed[LF_size//2:-LF_size//2],label='antenna 46')
pl1.plot(time[LF_size//2:-LF_size//2],amp_48_smoothed[LF_size//2:-LF_size//2],label='antenna 48')
pl1.plot(time[LF_size//2:-LF_size//2],amp_46_noise,label='noise 46')
pl1.plot(time[LF_size//2:-LF_size//2],amp_48_noise,label='noise 48')

pl1.plot(time[LF_size//2:-LF_size//2-noiseRange],amp_46_noise[:-noiseRange])
pl1.plot(time[LF_size//2:-LF_size//2-noiseRange],amp_48_noise[:-noiseRange])

pl1.set_ylim((-.5,1.))
pl1.set_title(dateTime + ', 4500 MHz, amp SNR %.3f, vis SNR %.3f'%(amp_46_smoothed[:-noiseRange].mean()/amp_46_noise[:-noiseRange].std(), vis_46_48_smoothed[:-noiseRange].mean()/vis_46_48_noise[:-noiseRange].std()))
pl1.set_xlabel('UTC')
pl1.set_ylabel('correlation')
pl1.grid()
pl1.legend()
PL.savefig(dateTime + '_4500_46_48.png')

PL.figure(figsize=(10,10))
pl1 = PL.subplot(111)
pl1.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl1.plot(time[LF_size//2:-LF_size//2],vis_101_104_smoothed[LF_size//2:-LF_size//2], label='visibility 101-104 amp')
pl1.plot(time[LF_size//2:-LF_size//2],vis_101_104_noise, label='visibility 101-104 noise')
pl1.plot(time[LF_size//2:-LF_size//2],amp_101_smoothed[LF_size//2:-LF_size//2],label='antenna 101')
pl1.plot(time[LF_size//2:-LF_size//2],amp_104_smoothed[LF_size//2:-LF_size//2],label='antenna 104')
pl1.plot(time[LF_size//2:-LF_size//2],amp_101_noise,label='noise 101')
pl1.plot(time[LF_size//2:-LF_size//2],amp_104_noise,label='noise 104')

pl1.plot(time[LF_size//2:-LF_size//2-noiseRange],amp_101_noise[:-noiseRange])
pl1.plot(time[LF_size//2:-LF_size//2-noiseRange],amp_104_noise[:-noiseRange])

pl1.set_ylim((-.5,1.))
pl1.set_title(dateTime + ', 4500 MHz, amp SNR %.3f, vis SNR %.3f'%(amp_101_smoothed[:-noiseRange].mean()/amp_101_noise[:-noiseRange].std(), vis_101_104_smoothed[:-noiseRange].mean()/vis_101_104_noise[:-noiseRange].std()))
pl1.set_xlabel('UTC')
pl1.set_ylabel('correlation')
pl1.grid()
pl1.legend()
PL.savefig(dateTime + '_4500_101_104.png')

