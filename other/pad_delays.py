#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 05:57:42 2019

@author: svlesovoi
"""
import pylab as PL
import numpy as NP

pads = [-37, 13, 15, 83, 2, 28, -16, -26, 31, 16, -6, -17, 47, 25, -22, -2,
        50, 45, 48, 44, 22, 66, 34, 39, -3, 32, 45, 32, 33, -28, -3, -36, 
        10, 20, 80, -16, -42, 60, -115, 30, 25, 24, -14, 89, 27, 32, 4, -25]
ants = ['192', '191', '190', '189', '188', '187', '186', '185', '184', '183', '182', '181', '180', '179', '178', '177',
        '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', 
        '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80']

pad_delays = dict(zip(ants, pads))

lAnts = list(pad_delays.keys())
lPads = list(pad_delays.values())

lfAnts = []
for s in range(len(lAnts)):
    lfAnts.append(int(lAnts[s]))

def ant_format(t, pos):
    if t >= 0. and t <= 47.:
        return  '%03d' % (lfAnts[int(t)])
    else: 
        return '%03d' % (t)

ind = NP.linspace(0,47,48, dtype='int')
alfAnts = NP.array(lAnts)

PL.figure()
pl0 = PL.subplot(111)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(ant_format))
pl0.xaxis.set_major_locator(PL.MultipleLocator(2))
pl0.set_xlim(0,47)
pl0.grid()

pl0.plot(ind, NP.array(lPads))

