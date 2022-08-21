#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 06:56:05 2020

@author: svlesovoi
"""
import numpy as np
import pylab as pl
from scipy import interpolate
x = np.arange(-5.01, 5.01, 0.1)
y = np.arange(-5.01, 5.01, 0.1)
xx, yy = np.meshgrid(x, y)
z = np.sin(xx**2+yy**2)
f = interpolate.interp2d(x, y, z, kind='cubic')

xnew = np.arange(-5.5, 5.5, .1)
ynew = np.arange(-5.5, 5.5, .1)
znew = f(xnew, ynew)
#pl.plot(x, z[:, 0], 'ro-', xnew, znew[:, 0], 'b-')
pl.show()
