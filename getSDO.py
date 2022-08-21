#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 09:00:02 2019

@author: sergey
"""

# Import the VSO client and create an instance
#from sunpy.net.vso import VSOClient
from sunpy.map import Map
import astropy.units as u
from sunpy.net import Fido, attrs as a

#client = VSOClient()
#query_response = client.query_legacy(tstart='2020/05/09 01:45:00', tend='2020/05/09 02:00:00', instrument='AIA', wave='304')
#query_response.sort(key=lambda x: x.time.start)
#print (query_response)
#results = client.get(query_response, site='rob')
#files = results.wait()

#results = Fido.search(a.Time("2020/07/18T00:00:00", "2020/07/18T07:00:00"), a.Instrument('AIA'),a.Wavelength(193*u.angstrom))
results = Fido.search(a.Time("2020/08/13T00:00:00", "2020/08/13T01:30:00"), a.Instrument('AIA'),a.Wavelength(193*u.angstrom))
files = Fido.fetch(results[0,1])

