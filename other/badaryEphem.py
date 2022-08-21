#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 08:27:47 2019

@author: sergey
"""

from astropy.coordinates import get_sun
import numpy as NP
from astropy import units as u
from astropy.time import Time
import pylab as PL
from astropy.coordinates import EarthLocation, AltAz
from optparse import OptionParser

def time2hhmm(t, pos):
    return T[0].iso.split(' ')[1].split('.')[0]
    return T[t].iso.split(' ')[1].split('.')[0]

def hm_format(t, pos):
  hh = int(t)
  t -= hh
  mm = int(t*60)
  return '%02d:%02d' % (hh,mm);

parser = OptionParser()
parser.add_option("-d", "--date", dest="givenDate", default = '2021-10-06')
parser.add_option("-s", "--start", dest="startTime", default = '00:00:00')
parser.add_option("-e", "--end", dest="endTime", default = '10:00:00')
parser.add_option("-n", "--size", dest="dTsize", default = 3600)
(badaryEphem_options, badaryEphem_args) = parser.parse_args()
givenDate = badaryEphem_options.givenDate
startTime = badaryEphem_options.startTime
endTime = badaryEphem_options.endTime
dTsize = int(badaryEphem_options.dTsize)

T_ = Time(givenDate)
T0 = Time(givenDate + ' ' + startTime)
T1 = Time(givenDate + ' ' + endTime)

dT = NP.linspace(0, (T1 - T0).sec / 3600, dTsize)*u.hour
dT_ = NP.linspace((T0 - T_).sec / 3600, (T1 - T_).sec / 3600, dTsize)*u.hour

T = T0 + dT

badary_rao = EarthLocation(lat=51.759371*u.deg, lon=102.219333*u.deg, height=799*u.m)
altazframe = AltAz(obstime=T, location=badary_rao)
sunaltazs = get_sun(T).transform_to(altazframe)

fig = PL.figure(figsize=(12,10))
fig.suptitle(givenDate)
pl = fig.subplots(nrows=1,ncols=1)
pl.xaxis.set_major_formatter(PL.FuncFormatter(hm_format))
#pl.xaxis.set_major_locator(PL.MultipleLocator(120))

pl.plot(dT_.to_value(), sunaltazs.az, label = 'azimuth')
pl.plot(dT_.to_value(), sunaltazs.alt, label = 'altitude')
pl.set_xlabel('time UTC')
pl.set_ylabel('degree')
pl.legend()
pl.grid()

saveEphemFile = open('badaryEphem_' + givenDate + '.txt','w')
saveEphemFile.write('date time altitude azimuth\n')
saveEphemFile.write('YYYY-MM-DD HH:MM:SS.MS deg deg\n')
for t in range(T.size):
    saveEphemFile.write('%s %.3f %.3f\n'%((T[t], sunaltazs.alt[t].to_value(), sunaltazs.az[t].to_value())))
    