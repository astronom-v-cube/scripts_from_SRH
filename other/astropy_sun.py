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
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

noon = Time('2021-09-30 05:01:06')
noon = Time('2021-10-17 04:56:0.')
#delta_noon = NP.linspace(-5, 5, 600)*u.hour
delta_noon = NP.linspace(-4.975, 4.975, 600)*u.hour
times = noon + delta_noon
badary_rao = EarthLocation(lat=51.759383*u.deg, lon=102.219309*u.deg, height=799*u.m)
#badary_rao = EarthLocation(lat=51.759371*u.deg, lon=102.219333*u.deg, height=799*u.m)
altazframe = AltAz(obstime=times, location=badary_rao)
sunaltazs = get_sun(times).transform_to(altazframe)

PL.figure()
PL.plot(delta_noon.to_value(), sunaltazs.az)
PL.plot(delta_noon.to_value(), sunaltazs.alt)
PL.grid()

# saveEphemFile = open('badaryEphem_20210930.txt','w')
# saveEphemFile.write('date time altitude azimuth\n')
# saveEphemFile.write('YYYY-MM-DD HH:MM:SS.MS deg deg\n')
# for t in range(delta_noon.size):
#     saveEphemFile.write('%s %.3f %.3f\n'%((times[t], sunaltazs.alt[t].to_value(), sunaltazs.az[t].to_value())))
    