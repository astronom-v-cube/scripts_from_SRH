#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 00:46:37 2020

@author: svlesovoi
"""

import pylab as PL
import numpy as NP
N = 512
L = 1000
alpha = NP.zeros(L)
rho = NP.zeros(L)

PL.figure()
for alpha2 in range(5):
    for i in range(L):
        alpha1 = (i * .001 + 1.)
        inputSignal = NP.random.randn(N)
        noise1 = NP.random.randn(N)
        noise2 = NP.random.randn(N)
        
        input1 = alpha1*inputSignal + noise1
        input2 = (alpha2 + 1.)*inputSignal + noise2
        alpha[i] = alpha1
        rho[i] = NP.mean(input1*input2)/(NP.std(input1)*NP.std(input2))
    PL.plot(alpha/(alpha2 + 1.),rho)
    

