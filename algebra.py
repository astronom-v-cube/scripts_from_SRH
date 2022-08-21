# -*- coding: utf-8 -*-
"""
Created on Sun Mar 18 02:16:48 2018

@author: Sergey
"""

import numpy as NP
import pylab as PL

x = NP.linspace(-10,10,1000)
y1 = - 2*x
y2 = 3 - x**2
PL.clf()
PL.plot(x,y1)
PL.plot(x,y2)
PL.grid()