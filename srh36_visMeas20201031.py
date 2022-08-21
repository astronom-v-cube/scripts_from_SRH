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

LF_size = 33
lf_filter = NP.ones(LF_size)

visScale = 1/(2e6*49)
ampScale = visScale / 128
#dateTime = '20201031T004345'
dateTime = '20201031T030324'
#dateTime = '20201031T041418'

amp_211 = fits.open('srh_amp_211_' + dateTime + '.fits')
amp_213 = fits.open('srh_amp_213_' + dateTime + '.fits')
amp_215 = fits.open('srh_amp_215_' + dateTime + '.fits')
amp_217 = fits.open('srh_amp_217_' + dateTime + '.fits')
amp_219 = fits.open('srh_amp_219_' + dateTime + '.fits')
amp_221 = fits.open('srh_amp_221_' + dateTime + '.fits')

vis_211_213 = fits.open('srh_vis_211_213_' + dateTime + '.fits')
vis_215_217 = fits.open('srh_vis_215_217_' + dateTime + '.fits')
vis_217_219 = fits.open('srh_vis_217_219_' + dateTime + '.fits')
vis_219_221 = fits.open('srh_vis_219_221_' + dateTime + '.fits')

time = vis_211_213[1].data['time']
global _timeStart
_timeStart = time[0]

vis_211_213_lcp = vis_211_213[1].data['vis_lcp']*visScale
vis_amp_211_213_lcp = NP.sqrt(NP.sin(NP.pi/2*vis_211_213_lcp.real)**2 + NP.sin(NP.pi/2*vis_211_213_lcp.imag)**2)
vis_211_213_rcp = vis_211_213[1].data['vis_rcp']*visScale
vis_amp_211_213_rcp = NP.sqrt(NP.sin(NP.pi/2*vis_211_213_rcp.real)**2 + NP.sin(NP.pi/2*vis_211_213_rcp.imag)**2)
vis_211_213_smoothed = NP.convolve((vis_amp_211_213_lcp + vis_amp_211_213_rcp)*.5, lf_filter, mode='same') / LF_size
vis_211_213_noise = vis_amp_211_213_lcp[LF_size//2:-LF_size//2] - vis_211_213_smoothed[LF_size//2:-LF_size//2]

vis_215_217_lcp = vis_215_217[1].data['vis_lcp']*visScale
vis_amp_215_217_lcp = NP.sqrt(NP.sin(NP.pi/2*vis_215_217_lcp.real)**2 + NP.sin(NP.pi/2*vis_215_217_lcp.imag)**2)
vis_215_217_rcp = vis_215_217[1].data['vis_rcp']*visScale
vis_amp_215_217_rcp = NP.sqrt(NP.sin(NP.pi/2*vis_215_217_rcp.real)**2 + NP.sin(NP.pi/2*vis_215_217_rcp.imag)**2)
vis_215_217_smoothed = NP.convolve((vis_amp_215_217_lcp + vis_amp_215_217_rcp)*.5, lf_filter, mode='same') / LF_size
vis_215_217_noise = vis_amp_215_217_lcp[LF_size//2:-LF_size//2] - vis_215_217_smoothed[LF_size//2:-LF_size//2]

vis_219_221_lcp = vis_219_221[1].data['vis_lcp']*visScale
vis_amp_219_221_lcp = NP.sqrt(NP.sin(NP.pi/2*vis_219_221_lcp.real)**2 + NP.sin(NP.pi/2*vis_219_221_lcp.imag)**2)
vis_219_221_rcp = vis_219_221[1].data['vis_rcp']*visScale
vis_amp_219_221_rcp = NP.sqrt(NP.sin(NP.pi/2*vis_219_221_rcp.real)**2 + NP.sin(NP.pi/2*vis_219_221_rcp.imag)**2)
vis_219_221_smoothed = NP.convolve((vis_amp_219_221_lcp + vis_amp_219_221_rcp)*.5, lf_filter, mode='same') / LF_size
vis_219_221_noise = vis_amp_219_221_lcp[LF_size//2:-LF_size//2] - vis_219_221_smoothed[LF_size//2:-LF_size//2]

vis_217_219_lcp = vis_217_219[1].data['vis_lcp']*visScale
vis_amp_217_219_lcp = NP.sqrt(NP.sin(NP.pi/2*vis_217_219_lcp.real)**2 + NP.sin(NP.pi/2*vis_217_219_lcp.imag)**2)
vis_217_219_rcp = vis_217_219[1].data['vis_rcp']*visScale
vis_amp_217_219_rcp = NP.sqrt(NP.sin(NP.pi/2*vis_217_219_rcp.real)**2 + NP.sin(NP.pi/2*vis_217_219_rcp.imag)**2)
vis_217_219_smoothed = NP.convolve((vis_amp_217_219_lcp + vis_amp_217_219_rcp)*.5, lf_filter, mode='same') / LF_size
vis_217_219_noise = vis_amp_217_219_lcp[LF_size//2:-LF_size//2] - vis_217_219_smoothed[LF_size//2:-LF_size//2]

time = time[LF_size//2:-LF_size//2]
vis_211_213_smoothed = vis_211_213_smoothed[LF_size//2:-LF_size//2]
vis_215_217_smoothed = vis_215_217_smoothed[LF_size//2:-LF_size//2]
vis_217_219_smoothed = vis_217_219_smoothed[LF_size//2:-LF_size//2]
vis_219_221_smoothed = vis_219_221_smoothed[LF_size//2:-LF_size//2]

noise_ind = [2150,2790]
noise_ind = [20150,20190]

SNR_211_213 = NP.abs(vis_211_213_smoothed[noise_ind[0]:noise_ind[1]].mean()/vis_211_213_noise[noise_ind[0]:noise_ind[1]].std())
SNR_215_217 = NP.abs(vis_215_217_smoothed[noise_ind[0]:noise_ind[1]].mean()/vis_215_217_noise[noise_ind[0]:noise_ind[1]].std())
SNR_217_219 = NP.abs(vis_217_219_smoothed[noise_ind[0]:noise_ind[1]].mean()/vis_217_219_noise[noise_ind[0]:noise_ind[1]].std())
SNR_219_221 = NP.abs(vis_219_221_smoothed[noise_ind[0]:noise_ind[1]].mean()/vis_219_221_noise[noise_ind[0]:noise_ind[1]].std())

freqList = amp_211[1].data['frequency']

PL.figure(figsize=(10,10))
pl0 = PL.subplot(111)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl0.plot(time,vis_211_213_smoothed,label='vis amp 211-213')
pl0.plot(time,vis_215_217_smoothed,label='vis amp 215-217')
pl0.plot(time,vis_217_219_smoothed,label='vis amp 217-219')
pl0.plot(time,vis_219_221_smoothed,label='vis amp 219-221')

pl0.set_title('pairs ' + dateTime + ', %d MHz' % freqList[0] + ' 219-221 interference!!!')
pl0.set_xlabel('UTC')
pl0.grid()
pl0.legend()

time = amp_211[1].data['time']
PL.figure(figsize=(10,10))
pl1 = PL.subplot(111)
pl1.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl1.plot(time,ampScale*.5*(amp_211[1].data['amp_lcp'] + amp_211[1].data['amp_rcp']),label='ch 08')
pl1.plot(time,ampScale*.5*(amp_213[1].data['amp_lcp'] + amp_213[1].data['amp_rcp']),label='ch 09')
pl1.plot(time,ampScale*.5*(amp_215[1].data['amp_lcp'] + amp_215[1].data['amp_rcp']),label='ch 10')
pl1.plot(time,ampScale*.5*(amp_217[1].data['amp_lcp'] + amp_217[1].data['amp_rcp']),label='ch 11')
pl1.plot(time,ampScale*.5*(amp_219[1].data['amp_lcp'] + amp_219[1].data['amp_rcp']),label='ch 12')
pl1.plot(time,ampScale*.5*(amp_221[1].data['amp_lcp'] + amp_221[1].data['amp_rcp']),label='ch 13')
pl1.set_title('single dish ' + dateTime + ', %d MHz' % freqList[0])
pl1.legend()

vis_vf_real = NP.sin(NP.pi/2*vis_215_217_lcp.real)
vis_vf_imag = NP.sin(NP.pi/2*vis_215_217_lcp.imag)
T = NP.linspace(0,vis_vf_imag.shape[0]-1,vis_vf_imag.shape[0])*2.8e-5 + 1.
