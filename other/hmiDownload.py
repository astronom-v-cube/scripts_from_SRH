#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 01:26:30 2021

@author: svlesovoi
"""

from sunpy.net import Fido, attrs as a
import sunpy.map

time = a.Time('2021/01/01', '2021/01/01')
series = a.jsoc.Series('hmi.synoptic_mr_polfil_720s')

result = Fido.search(time, series)
print(result)

files = Fido.fetch(result)
print(files)
