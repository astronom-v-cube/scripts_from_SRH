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

visScale = 1/(2e6*49)
ampScale = visScale / 128
LF_size = 19
lf_filter = NP.ones(LF_size)
#dateTime = '20201010T010027'
#dateTime = '20201010T020148'
dateTime = '20201010T042015'

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

amp_211_lcp = amp_211[1].data['amp_lcp'] * ampScale
amp_213_lcp = amp_213[1].data['amp_lcp'] * ampScale
vis_211_213_lcp = vis_211_213[1].data['vis_lcp'][0:sampleNumber] * visScale
vis_211_213_lcp_vf_real = NP.sin(NP.pi/2*vis_211_213_lcp.real)
vis_211_213_lcp_vf_imag = NP.sin(NP.pi/2*vis_211_213_lcp.imag)

amp_81_lcp = amp_81[1].data['amp_lcp'] * ampScale
amp_83_lcp = amp_83[1].data['amp_lcp'] * ampScale
vis_81_83_lcp = vis_81_83[1].data['vis_lcp'][0:sampleNumber] * visScale
vis_81_83_lcp_vf_real = NP.sin(NP.pi/2*vis_81_83_lcp.real)
vis_81_83_lcp_vf_imag = NP.sin(NP.pi/2*vis_81_83_lcp.imag)

amp_46_lcp = amp_46[1].data['amp_lcp'] * ampScale
amp_48_lcp = amp_48[1].data['amp_lcp'] * ampScale
vis_46_48_lcp = vis_46_48[1].data['vis_lcp'][0:sampleNumber] * visScale
vis_46_48_lcp_vf_real = NP.sin(NP.pi/2*vis_46_48_lcp.real)
vis_46_48_lcp_vf_imag = NP.sin(NP.pi/2*vis_46_48_lcp.imag)

amp_101_lcp = amp_101[1].data['amp_lcp'] * ampScale
amp_104_lcp = amp_104[1].data['amp_lcp'] * ampScale
vis_101_104_lcp = vis_101_104[1].data['vis_lcp'][0:sampleNumber] * visScale
vis_101_104_lcp_vf_real = NP.sin(NP.pi/2*vis_101_104_lcp.real)
vis_101_104_lcp_vf_imag = NP.sin(NP.pi/2*vis_101_104_lcp.imag)

PL.figure(figsize=(10,10))
pl0 = PL.subplot(111, )
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl0.plot(time,vis_211_213_lcp_vf_real,label='visibility 211-213 real')
pl0.plot(time,vis_211_213_lcp_vf_imag,label='visibility 211-213 imag')
pl0.plot(time,NP.sqrt(vis_211_213_lcp_vf_real**2 + vis_211_213_lcp_vf_imag**2), label='visibility 211-213 amp')
pl0.plot(time,amp_211_lcp,label='antenna 211')
pl0.plot(time,amp_213_lcp,label='antenna 213')
pl0.set_ylim((-.5,1.))
pl0.set_title(dateTime + ', 4500 MHz, real bias %.3f, imag bias %.3f'%(vis_211_213_lcp_vf_real.mean(), vis_211_213_lcp_vf_imag.mean()))
pl0.set_xlabel('UTC')
pl0.set_ylabel('correlation')
pl0.grid()
pl0.legend()
PL.savefig(dateTime + '_4500_211_213.png')

PL.figure(figsize=(10,10))
pl1 = PL.subplot(111)
pl1.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl1.plot(time,vis_81_83_lcp_vf_real,label='visibility 81-83 real')
pl1.plot(time,vis_81_83_lcp_vf_imag,label='visibility 81-83 imag')
pl1.plot(time,NP.sqrt(vis_81_83_lcp_vf_real**2 + vis_81_83_lcp_vf_imag**2), label='visibility 81-83 amp')
pl1.plot(time,amp_81_lcp,label='antenna 81')
pl1.plot(time,amp_83_lcp,label='antenna 83')
pl1.set_ylim((-.5,1.))
pl1.set_title(dateTime + ', 4500 MHz, real bias %.3f, imag bias %.3f'%(vis_81_83_lcp_vf_real.mean(), vis_81_83_lcp_vf_imag.mean()))
pl1.set_xlabel('UTC')
pl1.set_ylabel('correlation')
pl1.grid()
pl1.legend()
PL.savefig(dateTime + '_4500_81_83.png')

PL.figure(figsize=(10,10))
pl1 = PL.subplot(111)
pl1.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl1.plot(time,vis_46_48_lcp_vf_real,label='visibility 46-48 real')
pl1.plot(time,vis_46_48_lcp_vf_imag,label='visibility 46-48 imag')
pl1.plot(time,NP.sqrt(vis_46_48_lcp_vf_real**2 + vis_46_48_lcp_vf_imag**2), label='visibility 46-48 amp')
pl1.plot(time,amp_46_lcp,label='antenna 46')
pl1.plot(time,amp_48_lcp,label='antenna 48')
pl1.set_ylim((-.5,1.))
pl1.set_title(dateTime + ', 4500 MHz, real bias %.3f, imag bias %.3f'%(vis_46_48_lcp_vf_real.mean(), vis_46_48_lcp_vf_imag.mean()))
pl1.set_xlabel('UTC')
pl1.set_ylabel('correlation')
pl1.grid()
pl1.legend()
PL.savefig(dateTime + '_4500_46_48.png')

PL.figure(figsize=(10,10))
pl1 = PL.subplot(111)
pl1.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl1.plot(time,vis_101_104_lcp_vf_real,label='visibility 101-104 real')
pl1.plot(time,vis_101_104_lcp_vf_imag,label='visibility 101-104 imag')
pl1.plot(time,NP.sqrt(vis_101_104_lcp_vf_real**2 + vis_101_104_lcp_vf_imag**2), label='visibility 101-104 amp')
pl1.plot(time,amp_101_lcp,label='antenna 101')
pl1.plot(time,amp_104_lcp,label='antenna 104')
pl1.set_ylim((-.5,1.))
pl1.set_title(dateTime + ', 4500 MHz, real bias %.3f, imag bias %.3f'%(vis_101_104_lcp_vf_real.mean(), vis_101_104_lcp_vf_imag.mean()))
pl1.set_xlabel('UTC')
pl1.set_ylabel('correlation')
pl1.grid()
pl1.legend()
PL.savefig(dateTime + '_4500_101_104.png')

