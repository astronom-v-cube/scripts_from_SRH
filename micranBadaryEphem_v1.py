#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 10 08:27:47 2021

@author: sergey
"""

from astropy.coordinates import get_sun, get_moon
import numpy as NP
from astropy import units as u
from astropy.time import Time
import pylab as PL
from astropy.coordinates import EarthLocation, AltAz, SkyCoord

givenDate = '2022-01-02'
T0 = Time(givenDate + ' 00:00:00')
dT = 24 # in hours
dT = NP.linspace(0, dT, dT*3600)*u.hour

T = T0 + dT
badary_rao = EarthLocation(lat=51.759371*u.deg, lon=102.219333*u.deg, height=799*u.m)
altazframe = AltAz(obstime=T, location=badary_rao)
sunaltazs = get_sun(T).transform_to(altazframe)
moonaltazs = get_moon(T).transform_to(altazframe)
cas = SkyCoord.from_name('Cassiopeia A')
casaltazs = cas.transform_to(altazframe)

PL.figure()
PL.plot(dT.to_value(), sunaltazs.az, label = 'sun azimuth')
PL.plot(dT.to_value(), sunaltazs.alt, label = 'sun altitude')
PL.plot(dT.to_value(), moonaltazs.az, label = 'moon azimuth')
PL.plot(dT.to_value(), moonaltazs.alt, label = 'moon altitude')
PL.plot(dT.to_value(), casaltazs.az, label = 'Cas A azimuth')
PL.plot(dT.to_value(), casaltazs.alt, label = 'Cas A altitude')
PL.plot([dT[0].to_value(),dT[-1].to_value()], [180,180], label = 'meridium')
PL.xlabel('hour UTC')
PL.ylabel('degree')
PL.title(givenDate)
PL.legend()
PL.grid()

# saveEphemFile = open(givenDate +'_efmrd.txt','w')
# saveEphemFile.write('YYYY MM DD HH MM SS.sss DDD.dddddd DD.dddddd\n')
# for t in range(T.size):
#     reform_date = str(T[t]).replace(':',' ')
#     reform_date = reform_date.replace('-',' ')

#     # rotary azimuth on 180 degree
#     # if sunaltazs.az[t].to_value() - 180 > 0:
#     #     rotary_az = sunaltazs.az[t].to_value() - 180
#     # else:
#     #     rotary_az = 180 + sunaltazs.az[t].to_value() 
#     rotary_az = sunaltazs.az[t].to_value() + .25
#     if rotary_az - 180 > 0:
#         rotary_az = rotary_az - 180
#     else:
#         rotary_az = 180 + rotary_az

#     saveEphemFile.write('%s %.5f %.5f\n'%((reform_date, rotary_az, sunaltazs.alt[t].to_value())))
# saveEphemFile.close()
    