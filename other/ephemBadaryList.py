#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 02:06:04 2020

@author: svlesovoi
"""
import numpy as NP
from zeep import Client

def hhmmssFromSec(sec):
  hh = int(sec / 3600.)
  sec -= hh*3600.
  mm = int(sec / 60.)
  sec -= mm*60
  ss = int(sec)
  return '%02d:%02d:%02d' % (hh,mm,ss)

client = Client('http://ephemeris.rao.istp.ac.ru/?wsdl')

time = []
azimuth = []
altitude = []
ephemList = []

for seconds in NP.linspace(0,1*3600-60,1*60, dtype='int'):
    givenTime = hhmmssFromSec(seconds)
    result = client.service.Ephemeride('ssrt', 'sun', '2020-06-22T' + givenTime + 'Z')
    time.append(givenTime)
    azimuth.append(NP.deg2rad(float(result[0]['azimuth'])))
    altitude.append(NP.deg2rad(float(result[0]['altitude'])))
    ephemList.append('%s %f %f'%(givenTime, azimuth[-1], altitude[-1]))
    print(ephemList[-1])


