#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 08:52:49 2021

@author: sergey_lesovoi
"""

import numpy as NP
import pylab as PL
from astropy.coordinates import get_moon
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, solar_system_ephemeris

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

startDate = '2021-10-27'

solar_system_ephemeris.set('jpl')
start_t = Time(startDate + ' 00:00:00')
delta_t = NP.linspace(0, 10, 600)*u.hour
times = start_t + delta_t
badary_rao = EarthLocation(lat=51.759383*u.deg, lon=102.219309*u.deg, height=799*u.m)
altazframe = AltAz(obstime=times, location=badary_rao, pressure=0*u.hPa)#1013 hPa
sunaltazs = get_moon(times).transform_to(altazframe)


fig = PL.figure()
pl = fig.subplots(nrows=1,ncols=1)
fig.suptitle(startDate + ' coordinates of the Moon')
pl.set_xlabel('UTC')
pl.set_ylabel('degree')
pl.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
pl.xaxis.set_major_locator(PL.MultipleLocator(3600));
pl.plot(times.to_value(format='unix') - start_t.to_value(format='unix'),sunaltazs.az.to_value(),label='azimuth')
pl.plot(times.to_value(format='unix') - start_t.to_value(format='unix'),sunaltazs.alt.to_value(),label='altitude')
pl.grid()
pl.legend()

saveEphemFile = open('moon_ephem_badary_' + startDate + '.txt','w')
saveEphemFile.write('date time altitude azimuth\n')
saveEphemFile.write('YYYY-MM-DD HH:MM:SS.MS deg deg\n')
for t in range(delta_t.size):
    saveEphemFile.write('%s %.3f %.3f\n'%((times[t], sunaltazs.alt[t].to_value(), sunaltazs.az[t].to_value())))
