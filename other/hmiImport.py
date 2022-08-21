#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 01:12:32 2021

@author: svlesovoi
"""
from matplotlib import pyplot as plt
import sunpy.map
from sunpy.instr.aia import aiaprep
from sunpy.net import Fido, attrs as a

from astropy.coordinates import SkyCoord
from astropy import units as u

result = Fido.search(a.Time('2020-08-15T07:00:00', '2020-08-15T07:30:00'),
                     a.Instrument("aia"), a.Wavelength(171*u.angstrom),
                     a.vso.Sample(12*u.second))


file_download = Fido.fetch(result[0, 3], site='ROB')

aia1 = sunpy.map.Map(file_download[0])

aia = aiaprep(aia1)

plt.rc('font',family='serif')
plt.figure(figsize=[8,10])
aia.plot()
aia.draw_limb()
plt.grid(False)
plt.show()
