#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 01:36:25 2020

@author: svlesovoi
"""
import numpy as NP
import pylab as PL
from astropy.io import fits

visScale = 1/(2e6*49)
ampScale = visScale / 512
LF_size = 19
lf_filter = NP.ones(LF_size)

#vis_81_83 = fits.open('srh_vis_81_83_20200919T003309.fits')
#vis_211_213 = fits.open('srh_vis_211_213_20200919T003309.fits')
#vis_81_83 = fits.open('srh_vis_81_83_20200919T014601.fits')
#vis_211_213 = fits.open('srh_vis_211_213_20200919T014601.fits')
#vis_81_83 = fits.open('srh_vis_81_83_20200919T022401.fits')
#vis_211_213 = fits.open('srh_vis_211_213_20200919T022401.fits')
#vis_81_83 = fits.open('srh_vis_81_83_20200919T030545.fits')
#vis_211_213 = fits.open('srh_vis_211_213_20200919T030545.fits')
vis_81_83 = fits.open('srh_vis_81_83_20200919T052807.fits')
vis_211_213 = fits.open('srh_vis_211_213_20200919T052807.fits')

vis_81_83_lcp = vis_81_83[1].data['vis_lcp'] * visScale
vis_81_83_rcp = vis_81_83[1].data['vis_rcp'] * visScale
vis_81_83_lcp_vf_real = NP.sin(NP.pi/2*vis_81_83_lcp.real)
vis_81_83_rcp_vf_imag = NP.sin(NP.pi/2*vis_81_83_rcp.imag)

vis_211_213_lcp = vis_211_213[1].data['vis_lcp'] * visScale
vis_211_213_rcp = vis_211_213[1].data['vis_rcp'] * visScale
vis_211_213_lcp_vf_real = NP.sin(NP.pi/2*vis_211_213_lcp.real)
vis_211_213_rcp_vf_imag = NP.sin(NP.pi/2*vis_211_213_rcp.imag)

vis_81_83_real_smoothed = NP.convolve(vis_81_83_lcp.real, lf_filter, mode='same') / LF_size
vis_81_83_imag_smoothed = NP.convolve(vis_81_83_lcp.imag, lf_filter, mode='same') / LF_size
vis_81_83_abs_smoothed = NP.convolve(NP.abs(vis_81_83_lcp), lf_filter, mode='same') / LF_size
vis_211_213_real_smoothed = NP.convolve(vis_211_213_lcp.real, lf_filter, mode='same') / LF_size
vis_211_213_imag_smoothed = NP.convolve(vis_211_213_lcp.imag, lf_filter, mode='same') / LF_size
vis_211_213_abs_smoothed = NP.convolve(NP.abs(vis_211_213_lcp), lf_filter, mode='same') / LF_size

#amp_81 = fits.open('srh_amp_81_20200919T030545.fits')
#amp_211 = fits.open('srh_amp_211_20200919T030545.fits')
amp_81 = fits.open('srh_amp_81_20200919T052807.fits')
amp_83 = fits.open('srh_amp_83_20200919T052807.fits')
#amp_211 = fits.open('srh_amp_211_20200919T052807.fits')
amp_211 = fits.open('srh_amp_211_20201003T040018.fits')

amp_81_lcp = amp_81[1].data['amp_lcp'] * ampScale
amp_83_lcp = amp_83[1].data['amp_lcp'] * ampScale
amp_211_lcp = amp_211[1].data['amp_lcp'] * ampScale

amp_81_lcp_smoothed = NP.convolve(amp_81_lcp, lf_filter, mode='same') / LF_size
amp_83_lcp_smoothed = NP.convolve(amp_83_lcp, lf_filter, mode='same') / LF_size
amp_211_lcp_smoothed = NP.convolve(amp_211_lcp, lf_filter, mode='same') / LF_size

#freqList = vis_81_83[1].data['frequency']

PL.clf()
#PL.plot(NP.abs(vis_81_83_lcp))
#PL.plot(NP.abs(vis_81_83_lcp_vf))
#PL.plot(NP.abs(vis_211_213_lcp))
#PL.plot(NP.abs(vis_211_213_lcp_vf))
PL.plot(vis_81_83_real_smoothed)
PL.plot(vis_81_83_imag_smoothed)
PL.plot(NP.abs(vis_81_83_abs_smoothed))
PL.plot(vis_211_213_real_smoothed)
PL.plot(vis_211_213_imag_smoothed)
PL.plot(NP.abs(vis_211_213_abs_smoothed))
#PL.plot(vis_81_83_lcp.real - vis_81_83_smoothed)
#PL.plot(vis_211_213_lcp.real - vis_211_213_smoothed)
#PL.plot(amp_81_lcp_smoothed)
#PL.plot(amp_83_lcp_smoothed)
#PL.plot(amp_211_lcp_smoothed)
#PL.plot(amp_81_lcp - amp_81_lcp_smoothed)
#PL.plot(amp_211_lcp - amp_211_lcp_smoothed)
