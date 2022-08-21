#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 00:02:43 2019

@author: svlesovoi
"""

import numpy as NP
import pylab as PL

N = 1000

t = NP.linspace(0,1,N)
w = 2*NP.pi/.10

Z = NP.zeros(N, dtype='complex')
dZl = NP.zeros(N, dtype='complex')
dZr = NP.zeros(N, dtype='complex')

Z.real = NP.cos(w*t)
Z.imag = NP.sin(w*t)
dZl.real =  0.1*NP.cos(0.5*w*t)
dZl.imag =  0.1*NP.sin(0.5*w*t)
dZr.real =  0.1*NP.cos(0.5*w*t)
dZr.imag = -0.1*NP.sin(0.5*w*t)

idealPhase = NP.angle(Z)

Z *= (1 + 1.05*dZl.real)
#Z.imag = 1.1*Z.imag - 0.01

PL.clf()
PL.plot(t, Z.real)
PL.plot(t, Z.imag)
PL.plot(t, NP.abs(Z))
PL.plot(t, NP.angle(Z) - idealPhase)
