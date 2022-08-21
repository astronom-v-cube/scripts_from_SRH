#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 15:58:18 2019

@author: sergeylesovoi
"""

import numpy as NP
import pylab as PL
import scipy.signal

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

            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

CONS = 2

sAntPhaAccum = NP.zeros(16)
sAntPhaModel = NP.zeros(16)
for a in range(16):
    sAntPhaModel[a] = 180*(NP.random.rand() - .5)

sPha = NP.zeros(15)
for a in range(15):
    sPha[a] = 0 + sAntPhaModel[a] -  sAntPhaModel[a + 1] + 10*(NP.random.rand() - .5)

constraints = NP.zeros(CONS)
constraints[0] = 0#sAntPhaModel[0]
sAntPhaSolution, c, d, e = NP.linalg.lstsq(sPhaMatrix,NP.concatenate((sPha, constraints)))

sPhaLS = NP.zeros(15)
for a in range(15):
    sPhaLS[a] = 0 + sAntPhaSolution[a + 2] -  sAntPhaSolution[a + 1]

sAntPhaSolutionHand = NP.zeros(16)
sAntPhaSolutionHand[0] = 0#sAntPhaModel[0]
for ant in range(15):
    sAntPhaSolutionHand[ant + 1] = -sPha[ant] + sAntPhaSolutionHand[ant]

phaSlope = (sAntPhaModel[15] - sAntPhaModel[0])/15
phaSlopeArray = NP.zeros(16)
for ant in range(16):
    phaSlopeArray[ant] = ant*phaSlope

PL.clf()
PL.plot(sAntPhaModel, color='black', linewidth=3)
PL.plot(sAntPhaSolution[1:], 'o')
PL.plot(sAntPhaSolution[1:], color='green')
PL.plot(sAntPhaSolutionHand, 'o')
PL.plot(sAntPhaSolutionHand, color='red')

sAntPhaAccum = (NN*sAntPhaAccum + sAntPhaSolution[1:])/(NN+1)
NN += 1
PL.plot(sAntPhaAccum,color='magenta')
#
sigmaLS = NP.std(scipy.signal.detrend(sAntPhaSolution[1:] - sAntPhaModel))
sigmaAccum =  NP.std(scipy.signal.detrend(sAntPhaAccum - sAntPhaModel))
#
PL.title('diff std cur %f, diff std accum %f' % (sigmaLS, sigmaAccum))

