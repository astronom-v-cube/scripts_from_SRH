#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 18:19:44 2021

@author: sergeyvlesovoi
"""

import pylab as PL
import numpy as NP
import scipy.special
from scipy import signal

L=1024
arcsecRadius = 1020
degRadius = NP.deg2rad(arcsecRadius/3600)
u = NP.linspace(-L/2,L/2,1000)
v = NP.linspace(-L/2,L/2,1000)
j1_uArg = 2*NP.pi*degRadius*u
j1_vArg = 2*NP.pi*degRadius*v
xu,yv = NP.meshgrid(j1_uArg,j1_vArg)
apert = NP.zeros((1000,1000))
R = 60
for i in range(1000):
    for j in range(1000):
        if NP.sqrt((i-500)**2 + (j-500)**2) < R:
            apert[i,j] = 1

qSunSpectrum = NP.pi*scipy.special.jv(1,NP.sqrt(xu**2 + yv**2))/(NP.sqrt(xu**2 + yv**2))

win = signal.windows.triang(200)

uvPair = signal.fftconvolve(apert,apert,mode='same')/apert.sum()
pairQSunSpectrum = signal.fftconvolve(qSunSpectrum,uvPair,mode='same')/apert.sum()
PL.figure()
PL.plot(u,NP.abs(qSunSpectrum[500]))
PL.plot(u,NP.abs(pairQSunSpectrum[500]))
PL.grid()
