#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 15:58:18 2019

@author: sergeylesovoi
"""

import numpy as NP
import pylab as PL
import scipy.signal

sPhaMatrixSimple = [ \
            [1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1]]

sPhaMatrix = [ \
            [1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1], \

            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

CONS = 2

#sAntPhaModel = NP.zeros(16)
#for a in range(16):
#    sAntPhaModel[a] = 180*(NP.random.rand() - .5)

sPha = NP.zeros(15)
sunPha = -15.
constraints = NP.zeros(CONS)
constraints[0] = 0

K = 100
solutions = NP.zeros((16,K))
for k in range(K):
    for a in range(15):
        sPha[a] = sunPha + 20*(NP.random.rand() - .5) + sAntPhaModel[a] -  sAntPhaModel[a + 1]
    
    constraints[1] = NP.mean(sPha)
    sAntPhaSolution, c, d, e = NP.linalg.lstsq(sPhaMatrix,NP.concatenate((sPha, constraints)))
#    sAntPhaSolution, c, d, e = NP.linalg.lstsq(sPhaMatrixSimple,sPha)
    solutions[:,k] = sAntPhaSolution[1:]

PL.clf()
PL.plot(sAntPhaModel, color='black', linewidth=3)
#PL.plot(sAntPhaSolution[1:], 'o')
#PL.plot(sAntPhaSolution[1:], color='green')
PL.plot(solutions)

sigmaLS = NP.std(scipy.signal.detrend(sAntPhaSolution[1:] - sAntPhaModel))
PL.title('diff std cur %f' % (sigmaLS))

