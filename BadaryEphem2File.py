#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 10 08:27:47 2021

@author: sergey
"""

from astropy.coordinates import get_sun
import numpy as NP
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz
from optparse import OptionParser
import sys

parser = OptionParser()
parser.add_option("-d", "--date", dest="givenDate", default = '2021-10-06')
parser.add_option("-z", "--deltaAzimuth", dest="deltaAzimuth", default = '0')
parser.add_option("-a", "--deltaAltitude", dest="deltaAltitude", default = '0')
(badaryEphem_options, badaryEphem_args) = parser.parse_args()

if len(sys.argv) == 1:
    print('usage BadaryEphem2File.py -d 2021-11-14 -z 0.25 -a 0.\n or BadaryEphem2File.py --date 2021-11-14 --deltaAzimuth 0.25 -deltaAltitude 0')
else:    

    givenDate = badaryEphem_options.givenDate
    deltaAzimuth = float(badaryEphem_options.deltaAzimuth)
    deltaAltitude = float(badaryEphem_options.deltaAltitude)
    
    print('Calculate for %s, %f, %f'%(givenDate, deltaAzimuth, deltaAltitude))
    badary_rao = EarthLocation(lat=51.759371*u.deg, lon=102.219333*u.deg, height=799*u.m)
    T0 = Time(givenDate + ' 00:00:00',location=badary_rao)
    dT = 24 # in hours
    dT = NP.linspace(0, dT, dT*3600)*u.hour
    
    T = T0 + dT
    altazframe = AltAz(obstime=T, location=badary_rao)
    sunCoords = get_sun(T)
    sunaltazs = sunCoords.transform_to(altazframe)
    RA = sunCoords.ra
    LST = T.sidereal_time('mean')
    
    saveEphemFile = open(givenDate +'_efmrd.txt','w')
    saveEphemFile.write('YYYY MM DD HH MM SS.sss DDD.ddddd DD.ddddd\n')
    for t in range(T.size):
        reform_date = str(T[t]).replace(':',' ')
        reform_date = reform_date.replace('-',' ')
    
        cur_az  = sunaltazs.az[t].to_value() + deltaAzimuth
        cur_alt = sunaltazs.alt[t].to_value() + deltaAltitude
        if cur_az - 180 > 0:
            cur_az = cur_az - 180
        else:
            cur_az = 180 + cur_az
    
        curHourAngle = LST[t].to_value()*15-RA[t].to_value()
        saveEphemFile.write('%s %.5f %.5f %.5f\n'%((reform_date, cur_az, cur_alt, curHourAngle)))
    saveEphemFile.close()
