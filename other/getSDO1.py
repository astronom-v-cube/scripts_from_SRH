#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 10:26:04 2019

@author: maria
"""

import astropy.units as u
from sunpy.net import Fido, attrs as a

attrs_time = a.Time('2022/02/15 03:30', '2022/02/15 03:35')
result = Fido.search(attrs_time, a.Instrument('aia'), a.Wavelength(304*u.angstrom))

downloaded_files = Fido.fetch(result)
print(downloaded_files)

result = Fido.search(attrs_time, a.Instrument('aia'), a.Wavelength(171*u.angstrom))

downloaded_files = Fido.fetch(result)
print(downloaded_files)

result = Fido.search(attrs_time, a.Instrument('hmi'), a.Physobs.los_magnetic_field)

downloaded_files = Fido.fetch(result)
print(downloaded_files)
