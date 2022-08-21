#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 10:00:56 2021

@author: sergeyvlesovoi
"""

from astropy.io import fits
import numpy as NP
import pylab as PL

def quiteSunSfuAsFrequncy(freq, freq0):
    return 80 + ((freq - freq0)*2E-6)**2

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

#times = [13600,14200]
times = [0,20000]
sky_times = [1000,1500]
sun_times = [26400,26800]
fluxAndCp = fits.open('srh_cp_20210703.fits')
frequencies = fluxAndCp[1].data['frequencies']

#flareCp = fluxAndCp[2].data['I'][:,times[0]:times[1]].copy()
#flareFlux = fluxAndCp[2].data['flux_I'][:,times[0]:times[1]].copy()
flareTime = fluxAndCp[2].data['time']
flareCp = fluxAndCp[2].data['I'][:].copy()
flareFlux = fluxAndCp[2].data['flux_I'][:].copy()

meanSky = []
meanSun = []

for f in range(frequencies.shape[0]):
    sunMean = fluxAndCp[2].data['flux_I'][f,sun_times].mean()
    skyMean = sunMean/2
    meanSky.append(skyMean)
    meanSun.append(sunMean)

for f in range(frequencies.shape[0]):
    flareFlux[f,:] -= meanSky[f]
    flareFlux[f,:] /= (meanSun[f] - meanSky[f])
    flareFlux[f,:] *= quiteSunSfuAsFrequncy(frequencies[f], frequencies[0]) 

meanSky = []
meanSun = []

for f in range(frequencies.shape[0]):
    meanSky.append(fluxAndCp[2].data['I'][f,sky_times].mean())
    meanSun.append(fluxAndCp[2].data['I'][f,sun_times].mean())

#for f in range(frequencies.shape[0]):
#    flareCp[f,:] -= meanSky[f]
#    flareCp[f,:] /= (meanSun[f] - meanSky[f])
#    flareCp[f,:] *= quiteSunSfuAsFrequncy(frequencies[f], frequencies[0]) 
#
#for f in range(frequencies.shape[0]):
#    flareFlux[f,:] -= flareFlux[f,0]
#    flareCp[f,:] -= flareCp[f,0]

#flareTb = flareCp + flareFlux
#flareSize = flareFlux/flareTb
#PL.figure()
#PL.xlim(flareTime[0,19100],flareTime[0,19700])
#for f in range(frequencies.shape[0]):
#    PL.plot(flareTime[0],(flareCp[0]-flareCp[0,19330])*0.7e9+Tbmax_0703_2800[0],'.',color='red',markersize=0.3)
#    PL.plot(flareTime[1],(flareCp[1]-flareCp[1,19330])*0.7e9+Tbmax_0703_3100[0],'.',color='green',markersize=0.3)
#    PL.plot(flareTime[2],(flareCp[2]-flareCp[2,19330])*0.7e9+Tbmax_0703_3400[0],'.',color='blue',markersize=0.3)

PL.figure()
PL.title('flux')
PL.xlim(flareTime[0,19100],flareTime[0,19700])
PL.plot(flareTime[0],flareFlux[0]-flareFlux[0,19330]*0,color='red',markersize=0.3)
PL.plot(flareTime[1],flareFlux[1]-flareFlux[1,19330]*0,color='green',markersize=0.3)
PL.plot(flareTime[2],flareFlux[2]-flareFlux[2,19330]*0,color='blue',markersize=0.3)
PL.plot(time_0703_2800,(flux_0703_2800 - flux_0703_2800[0]*0)*1.,'.')
PL.plot(time_0703_3100,(flux_0703_3100 - flux_0703_3100[0]*0)*1.,'.')
PL.plot(time_0703_3400,(flux_0703_3400 - flux_0703_3400[0]*0)*1.,'.')
PL.plot(time_0716_2800,(flux_0716_2800 - flux_0703_2800[0]*0)*1.,'.')
PL.plot(time_0716_3100,(flux_0716_3100 - flux_0703_3100[0]*0)*1.,'.')
PL.plot(time_0716_3400,(flux_0716_3400 - flux_0703_3400[0]*0)*1.,'.')
PL.plot(time_0718_2800,(flux_0718_2800 - flux_0703_2800[0]*0)*1.,'.')
PL.plot(time_0718_3100,(flux_0718_3100 - flux_0703_3100[0]*0)*1.,'.')
PL.plot(time_0718_3400,(flux_0718_3400 - flux_0703_3400[0]*0)*1.,'.')

commonFlux = (flareFlux[0]+flareFlux[1]+flareFlux[2])/3.
PL.rcParams.update({'font.size': 14})
fig = PL.figure()
fig.suptitle('SRH 20210703, 2.8-3.4 GHz')
sub = fig.add_subplot(1,1,1);
sub.set_ylabel('s.f.u.');
sub.set_xlabel('UT');
sub.xaxis.set_major_locator(PL.MultipleLocator(300));
sub.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
sub.xaxis.set_minor_locator(PL.MultipleLocator(20));
sub.set_xlim(flareTime[0,19100],flareTime[0,19700])
sub.set_ylim(0,100)
sub.grid()
sub.plot(flareTime[0],commonFlux - commonFlux[19100],color='black',markersize=0.3,label='from single antenna')
sub.plot(time_0703_2800,(flux_0703_2800 - flux_0703_2800[0]*2.0 + flux_0703_3100 + flux_0703_3400)/3.*1.75,'.',color='green',label='from images')
sub.plot(time_0716_2800,(flux_0716_2800 - flux_0703_2800[0]*0.7 + flux_0716_3100 + flux_0716_3400)/3.*1.3,'.',color='green')
sub.plot(time_0718_2800,(flux_0718_2800 - flux_0703_2800[0]*0.7 + flux_0718_3100 + flux_0718_3400)/3.*1.15,'.',color='green')
sub.legend()

fig = PL.figure()
fig.suptitle('SRH 20210703, 2.8-3.4 GHz')
sub = fig.add_subplot(1,1,1);
sub.set_ylabel('brightness temperature');
sub.set_xlabel('UT');
sub.xaxis.set_major_locator(PL.MultipleLocator(300));
sub.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
sub.xaxis.set_minor_locator(PL.MultipleLocator(20));
sub.set_xlim(flareTime[0,19100],flareTime[0,19700])
sub.set_ylim(0,2.5e7)
sub.grid()

Tbmax_0703_2800 = NP.array(Tbmax_0703_2800)
Tbmax_0703_3100 = NP.array(Tbmax_0703_3100)
Tbmax_0703_3400 = NP.array(Tbmax_0703_3400)
Tbmax_0703 = (Tbmax_0703_2800 + Tbmax_0703_3100 + Tbmax_0703_3400)/3

Tbmax_0716_2800 = NP.array(Tbmax_0716_2800)
Tbmax_0716_3100 = NP.array(Tbmax_0716_3100)
Tbmax_0716_3400 = NP.array(Tbmax_0716_3400)
Tbmax_0716 = (Tbmax_0716_2800 + Tbmax_0716_3100 + Tbmax_0716_3400)/3

Tbmax_0718_2800 = NP.array(Tbmax_0718_2800)
Tbmax_0718_3100 = NP.array(Tbmax_0718_3100)
Tbmax_0718_3400 = NP.array(Tbmax_0718_3400)
Tbmax_0718 = (Tbmax_0718_2800 + Tbmax_0718_3100 + Tbmax_0718_3400)/3

sub.plot(flareTime[0],(flareCp[0]-flareCp[0,19330])*0.8e9 + Tbmax_0703_2800[0],color='red',label='from correaltion plot')
sub.plot(time_0703_2800,Tbmax_0703,'.',color='green',label='from images')
sub.plot(time_0716_2800,Tbmax_0716,'.',color='green')
sub.plot(time_0718_2800,Tbmax_0718,'.',color='green')
sub.legend()
