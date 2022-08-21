#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 10 08:27:47 2021

@author: sergey
"""

from astropy.coordinates import get_moon
import numpy as NP
from astropy import units as u
from astropy.time import Time
import pylab as PL
from astropy.coordinates import EarthLocation, AltAz

givenDate = '2022-03-23'
T0 = Time(givenDate + ' 00:00:00')
dT = 24 # in hours
dT = NP.linspace(0, dT, dT*3600)*u.hour

T = T0 + dT
badary_rao = EarthLocation(lat=51.759371*u.deg, lon=102.219333*u.deg, height=799*u.m)
altazframe = AltAz(obstime=T, location=badary_rao)
sunaltazs = get_moon(T).transform_to(altazframe)

PL.figure()
PL.plot(dT.to_value(), sunaltazs.az, label = 'azimuth')
PL.plot(dT.to_value(), sunaltazs.alt, label = 'altitude')
PL.xlabel('hour UTC')
PL.ylabel('degree')
PL.title(givenDate)
PL.legend()
PL.grid()

saveEphemFile = open(givenDate +'_efmrd_moon.txt','w')
saveEphemFile.write('YYYY MM DD HH MM SS.sss DDD.dddddd DD.dddddd\n')
for t in range(T.size):
    reform_date = str(T[t]).replace(':',' ')
    reform_date = reform_date.replace('-',' ')

    # rotary azimuth on 180 degree
    if sunaltazs.az[t].to_value() - 180 > 0:
        rotary_az = sunaltazs.az[t].to_value() - 180
    else:
        rotary_az = 180 + sunaltazs.az[t].to_value() 

    saveEphemFile.write('%s %.5f %.5f\n'%((reform_date, rotary_az, sunaltazs.alt[t].to_value())))
saveEphemFile.close()
    