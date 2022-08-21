#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 08:44:00 2021

@author: sergey_lesovoi
"""
import numpy as NP
import pylab as PL
from astropy.coordinates import get_sun
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, solar_system_ephemeris

def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]

def hhmm2s(hhmm):
    text = hhmm.split(':')
    return int(text[0])*3600 + int(text[1])*60

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

startDate = '2021-10-17'
#------------------------------------------------------------------------------
sun_eph_JPL_file = open('sun_eph_JPL.txt')
sun_eph_JPL_text = sun_eph_JPL_file.readlines()

times_JPL = []
az_JPL = []
alt_JPL = []

for i in NP.arange(4,604):
    fields = sun_eph_JPL_text[i].split(' ')
    fields = remove_values_from_list(fields,'')
    if len(fields) == 11:
        times_JPL.append(hhmm2s(fields[1]))
        az_JPL.append(float(fields[9]))
        alt_JPL.append(float(fields[10].split('\n')[0]))

times_JPL = NP.array(times_JPL)
az_JPL = NP.array(az_JPL)
alt_JPL = NP.array(alt_JPL)

#------------------------------------------------------------------------------
sun_eph_IAA_file = open('sun_eph_IAA.txt')
sun_eph_IAA_text = sun_eph_IAA_file.readlines()

times_IAA = []
az_IAA = []
alt_IAA = []

for i in NP.arange(13,613):
    fields = sun_eph_IAA_text[i].split(' ')
    fields = remove_values_from_list(fields,'')
    times_IAA.append(hhmm2s(fields[3] + ':' + fields[4] + ':' + fields[5]))
    az_IAA.append(float(fields[6]))
    if (az_IAA[-1] > 180):
        az_IAA[-1] -= 180
    else:
        az_IAA[-1] += 180
    alt_IAA.append(float(fields[7].split('\n')[0]))

times_IAA = NP.array(times_IAA)
az_IAA = NP.array(az_IAA)
alt_IAA = NP.array(alt_IAA)

#------------------------------------------------------------------------------
solar_system_ephemeris.set('jpl')
start_t = Time(startDate + ' 00:00:00')
delta_t = NP.linspace(0, 10, 600)*u.hour
times = start_t + delta_t
badary_rao = EarthLocation(lat=51.759383*u.deg, lon=102.219309*u.deg, height=799*u.m)
altazframe = AltAz(obstime=times, location=badary_rao, pressure=0*u.hPa)#1013 hPa
sunaltazs = get_sun(times).transform_to(altazframe)
 
#------------------------------------------------------------------------------
fig = PL.figure()
pl = fig.subplots(nrows=1,ncols=1)
fig.suptitle(startDate + '\nDifferences beteween IAA, astropy (ERFA) and Horizons (JPL) Sun ephemeris')
pl.set_xlabel('UTC')
pl.set_ylabel('degree')
pl.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl.xaxis.set_major_locator(PL.MultipleLocator(3600));
# PL.plot(times_JPL,az_JPL)
# PL.plot(times_JPL,alt_JPL)

#PL.plot(times_IAA, sunaltazs.az.to_value())
#PL.plot(times_IAA, sunaltazs.alt.to_value())

#PL.plot(times_IAA,az_IAA)
#PL.plot(times_IAA,alt_IAA)

#PL.plot(times_IAA, az_IAA - az_JPL,label = 'Az IAA - JPL')
#PL.plot(times_IAA, alt_IAA - alt_JPL,label = 'Alt IAA - JPL')

pl.plot(times_IAA, az_IAA - az_JPL,label='Az IAA-JPL')
pl.plot(times_IAA, alt_IAA - alt_JPL, label = 'Alt IAA-JPL')

pl.plot(times_IAA, sunaltazs.az.to_value() - az_JPL,label='Az ERFA-JPL')
pl.plot(times_IAA, sunaltazs.alt.to_value() - alt_JPL, label = 'Alt ERFA-JPL')

pl.legend()
pl.grid()


