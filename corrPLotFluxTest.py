#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 13:05:47 2021

@author: sergeyvlesovoi
"""


import numpy as NP
import pylab as PL

N = 1000
M = 20
L = N // 5
x = NP.zeros(N)
Tb = []
theta = []
flux = []
corr = []
for i in range(M):
    x[:] = 0
    x[N//2 - i:N//2 + i] = 1
    y = NP.abs(NP.fft.fft(x))
    Tb.append(1)
    theta.append(i)
    flux.append(y[0])
    corr.append(y[1:L].mean())

PL.figure()    
PL.plot(flux)
PL.plot(corr)

Tb = []
theta = []
flux = []
corr = []
for i in range(M):
    x[:] = 0
    x[N//2 - (M-i):N//2 + (M-i)] = NP.exp(0.2*i)
    y = NP.abs(NP.fft.fft(x))
    Tb.append(NP.exp(0.2*i))
    theta.append(M-i)
    flux.append(y[0])
    corr.append(y[1:L].mean())

PL.figure()    
PL.plot(flux)
PL.plot(corr)