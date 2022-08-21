#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 02:58:22 2020

@author: svlesovoi
"""

import numpy as NP
import pylab as PL
from zeep import Client

def hhmmss_format(t):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  t -= mm*60.;
  ss = (int)(t);
  return '%02d:%02d:%02d' % (hh,mm,ss);

client = Client('http://ephemeris.rao.istp.ac.ru/?wsdl')

azimuth = []
altitude = []
selected_t = []
selected_azimuth = []
selected_altitude = []

t0 = 2*3600
N = 360
given_t = NP.linspace(t0,t0 + (N - 1)*10,N) 

for t in given_t:
    dateTimeText = hhmmss_format(t)
    result = client.service.Ephemeride('ssrt', 'sun', '2020-06-22T' + dateTimeText + 'Z')
    azimuth.append(float(result[0]['azimuth']))
    altitude.append(float(result[0]['altitude']))
    if (t % 300 == 0. or t == t0 + (N - 1)*10):
        selected_t.append(t)
        selected_azimuth.append(azimuth[-1])
        selected_altitude.append(altitude[-1])

int_azimuth = NP.interp(given_t, selected_t, selected_azimuth)
int_altitude = NP.interp(given_t, selected_t, selected_altitude)

PL.clf()

PL.plot(given_t, azimuth)
PL.plot(selected_t, selected_azimuth, '.')
PL.plot(given_t, int_azimuth)

PL.plot(given_t, altitude)
PL.plot(selected_t, selected_altitude, '.')
PL.plot(given_t, int_altitude)

print(NP.std((azimuth - int_azimuth)*60))
print(NP.std((altitude - int_altitude)*60))



