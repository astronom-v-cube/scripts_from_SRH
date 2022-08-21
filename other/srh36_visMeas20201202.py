#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 01:36:25 2020

@author: svlesovoi
"""
import numpy as NP
import pylab as PL
from astropy.io import fits
import scipy.signal

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
dateTime = '20201202T024031'

vis_211_213 = fits.open('srh_vis_211_213_' + dateTime + '.fits')
vis_213_215 = fits.open('srh_vis_213_215_' + dateTime + '.fits')
vis_215_217 = fits.open('srh_vis_215_217_' + dateTime + '.fits')
vis_217_219 = fits.open('srh_vis_217_219_' + dateTime + '.fits')

vis_211_215 = fits.open('srh_vis_211_215_' + dateTime + '.fits')
vis_213_217 = fits.open('srh_vis_213_217_' + dateTime + '.fits')
vis_215_219 = fits.open('srh_vis_215_219_' + dateTime + '.fits')

time = vis_211_213[1].data['time']
global _timeStart
_timeStart = time[0]

vis_211_213_lcp = vis_211_213[1].data['vis_lcp']*visScale
vis_211_213_rcp = vis_211_213[1].data['vis_rcp']*visScale

vis_213_215_lcp = vis_213_215[1].data['vis_lcp']*visScale
vis_213_215_rcp = vis_213_215[1].data['vis_rcp']*visScale

vis_215_217_lcp = vis_215_217[1].data['vis_lcp']*visScale
vis_215_217_rcp = vis_215_217[1].data['vis_rcp']*visScale

vis_217_219_lcp = vis_217_219[1].data['vis_lcp']*visScale
vis_217_219_rcp = vis_217_219[1].data['vis_rcp']*visScale

vis_211_215_lcp = vis_211_215[1].data['vis_lcp']*visScale
vis_211_215_rcp = vis_211_215[1].data['vis_rcp']*visScale

vis_213_217_lcp = vis_213_217[1].data['vis_lcp']*visScale
vis_213_217_rcp = vis_213_217[1].data['vis_rcp']*visScale

vis_215_219_lcp = vis_215_219[1].data['vis_lcp']*visScale
vis_215_219_rcp = vis_215_219[1].data['vis_rcp']*visScale

freqList = vis_211_213[1].data['frequency']

pha_211_213 = (NP.unwrap(NP.angle(vis_211_213_lcp[200:],deg=False)))
pha_213_215 = (NP.unwrap(NP.angle(vis_213_215_lcp[200:],deg=False)))
pha_215_217 = (NP.unwrap(NP.angle(vis_215_217_lcp[200:],deg=False)))
pha_217_219 = (NP.unwrap(NP.angle(vis_217_219_lcp[200:],deg=False)))

pha_211_215 = (NP.unwrap(NP.angle(vis_211_215_lcp[200:],deg=False)))
pha_213_217 = (NP.unwrap(NP.angle(vis_213_217_lcp[200:],deg=False)))
pha_215_219 = (NP.unwrap(NP.angle(vis_215_219_lcp[200:],deg=False)))

piScale1 = 1
piScale2 = 1
piScale3 = 1

PL.figure(figsize=(10,10))
pl0 = PL.subplot(111)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl0.plot(time[200:],NP.rad2deg(pha_211_213 + pha_213_215 - pha_211_215 + 2*NP.pi*piScale1),label='211,213,215')
pl0.plot(time[200:],NP.rad2deg(pha_213_215 + pha_215_217 - pha_213_217 + 2*NP.pi*piScale2),label='213,215,217')
pl0.plot(time[200:],NP.rad2deg(pha_215_217 + pha_217_219 - pha_215_219 + 2*NP.pi*piScale3),label='215,217,219')
pl0.plot(time[200:],NP.rad2deg(pha_211_213)*.05,label='211,213 *.05')
pl0.plot(time[200:],NP.rad2deg(pha_213_215)*.05,label='213,215 *.05')
pl0.plot(time[200:],NP.rad2deg(pha_211_215)*.05,label='211,215 *.05')

pl0.set_title('closure phase ' + dateTime + ', %d MHz' % freqList[0])
pl0.set_ylim([-180,180])
pl0.set_xlabel('UTC')
pl0.set_ylabel('degree')
pl0.grid()
pl0.legend()

#PL.figure(figsize=(10,10))
#pl1 = PL.subplot(111)
#pl1.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
#pl1.plot(time[200:],NP.rad2deg(pha_211_213),label='211,213')
#pl1.plot(time[200:],NP.rad2deg(pha_213_215),label='213,215')
#pl1.plot(time[200:],NP.rad2deg(pha_211_215),label='211,215')
#
#pl1.set_title('phase ' + dateTime + ', %d MHz' % freqList[0])
#pl1.set_ylim([-18000,18000])
#pl1.set_xlabel('UTC')
#pl1.set_ylabel('degree')
#pl1.grid()
#pl1.legend()

time = vis_211_213[1].data['time']
vis_vf_real = NP.sin(NP.pi/2*vis_215_217_lcp.real)
vis_vf_imag = NP.sin(NP.pi/2*vis_215_217_lcp.imag)
T = NP.linspace(0,vis_vf_imag.shape[0]-1,vis_vf_imag.shape[0])*2.8e-5 + 1.
