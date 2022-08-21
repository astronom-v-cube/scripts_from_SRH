#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 02:06:04 2020

@author: svlesovoi
"""
import numpy as NP
from zeep import Client

client = Client('http://ephemeris.rao.istp.ac.ru/?wsdl')
result = client.service.Ephemeride('ssrt', 'sun', '2020-06-22T05:14:00Z')
azimuth = NP.deg2rad(float(result[0]['azimuth']))
altitude = NP.deg2rad(float(result[0]['altitude']))
declination = NP.deg2rad(float(result[0]['declination']))
sinHA = NP.sin(azimuth)*NP.cos(altitude)/NP.cos(declination)
print (result, NP.arcsin(sinHA))

ephList = client.service.ListEph('2020-06-22T05:14:00Z', '2020-06-22T06:14:00Z')
ephList = ephList[0]['DT'].replace('),','').split('(')[1:]
dateAltAzList = []
for i in range(len(ephList)):
    dateAltAzList.append(ephList[i].replace('\'','').split(','))


