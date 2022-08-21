#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 08:23:06 2022

@author: sergeyvlesovoi
"""


import pylab as PL
import numpy as NP

def quantizer(inputSignal, q = 1, levels = 1,quantizerType = 'even'):
    if quantizerType == 'even':
        return NP.clip(inputSignal//q*q + q/2,-(levels - .5)*q,(levels - .5)*q)
    elif quantizerType == 'evenZero':
        res = NP.clip(inputSignal//q*q + q/2,-(levels - .5)*q,(levels - .5)*q)
        res[NP.where(inputSignal == 0)[0]] = 0
        return res
    else:
        return NP.clip((inputSignal + q/2)//q*q,-levels*q,levels*q)

    
L = 500
N = 1000
M = 1
scaling = 1
q = 500
levels = 7

alphas = NP.linspace(0,1.,L)

autoCorr = NP.zeros((M,L))
autoCorrADC = NP.zeros((M,L))
autoCorrOdd = NP.zeros((M,L))
autoCorrEven = NP.zeros((M,L))
autoCorrEvenZero = NP.zeros((M,L))

zeroCountADC = NP.zeros((M,L))
zeroCountOdd = NP.zeros((M,L))
zeroCountEven = NP.zeros((M,L))
zeroCountEvenZero = NP.zeros((M,L))

for i in range(M):
    for j in range(L):
        alpha = alphas[j]
        inputSignal = (alpha*NP.random.randn(N)*scaling).astype(float)
        inputSignalADC = quantizer(inputSignal,q=1e-3,quantizerType='odd',levels=2047)
        inputSignalOdd = quantizer(inputSignalADC,q=1,quantizerType='odd',levels=7)
        inputSignalEven = quantizer(inputSignalADC,q=1,quantizerType='even',levels=7)
        inputSignalEvenZero = quantizer(inputSignalADC,q=1,quantizerType='evenZero',levels=7)
        
        zeroCountADC[i,j] = len(NP.where(inputSignalADC == 0)[0])
        zeroCountOdd[i,j] = len(NP.where(inputSignalOdd == 0)[0])
        zeroCountEven[i,j] = len(NP.where(inputSignalEven == 0)[0])
        zeroCountEvenZero[i,j] = len(NP.where(inputSignalEvenZero == 0)[0])

        autoCorr[i,j] = NP.mean(inputSignal*inputSignal)
        autoCorrADC[i,j] = NP.mean(inputSignalADC*inputSignalADC)
        autoCorrOdd[i,j] = NP.mean(inputSignalOdd*inputSignalOdd)
        autoCorrEven[i,j] = NP.mean(inputSignalEven*inputSignalEven)
        autoCorrEvenZero[i,j] = NP.mean(inputSignalEvenZero*inputSignalEvenZero)
    
fig = PL.figure()
fig.suptitle(r'$\sigma = %d, $q = %d'%(scaling,q))
pl = fig.subplots(nrows = 1, ncols = 5)
pl[0].plot(alphas,autoCorr[0],'.',label='float')
pl[1].plot(alphas,autoCorrADC[0],'.',label='ADC')
pl[2].plot(alphas,autoCorrOdd[0],'.',label='odd')
pl[3].plot(alphas,autoCorrEven[0],'.',label='even')
pl[4].plot(alphas,autoCorrEvenZero[0],'.',label='evenZero')

for i in range(5):
    pl[i].set_ylim(0,1.5)
    pl[i].set_xlabel('input sigma')
    pl[i].set_ylabel('output power')
    pl[i].grid()
    pl[i].legend()

PL.figure()
PL.plot(alphas,zeroCountADC[0]/N*100,label='ADC')
PL.plot(alphas,zeroCountOdd[0]/N*100,label='odd')
PL.plot(alphas,zeroCountEven[0]/N*100,label='even')
PL.plot(alphas,zeroCountEvenZero[0]/N*100,label='evenZero')
PL.xlabel('input sigma')
PL.ylabel('zero count')
PL.legend()
PL.grid()