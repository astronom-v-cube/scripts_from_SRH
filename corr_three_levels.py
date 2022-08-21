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
    else:
        return NP.clip((inputSignal + q/2)//q*q,-levels*q,levels*q)

    
N = 500
L = 500

alphas = NP.linspace(0,10,L)
rhos = NP.zeros(L)
threeRhos = NP.zeros(L)
covs = NP.zeros(L)
threeCovs = NP.zeros(L)
autocorrs = NP.zeros(L)
oddThreeAutocorrs = NP.zeros(L)
evenThreeAutocorrs = NP.zeros(L)
oddNNcovs = NP.zeros(L)
NNcovs = NP.zeros(L)

scaling = 1000
q = 100
levels = 7
for j in range(L):
    alpha = alphas[j]
    inputSignal = (NP.random.randn(N)*scaling).astype(int)
    noise1 = (NP.random.randn(N)*scaling).astype(int)
    noise2 = (NP.random.randn(N)*scaling).astype(int)
    input1 = (NP.sqrt(alpha)*inputSignal + noise1).astype(int)
    input2 = (NP.sqrt(alpha)*inputSignal + noise2).astype(int)
    oddThreeLevelsInput1 = quantizer(input1,q=q,quantizerType='odd',levels=levels)
    oddThreeLevelsInput2 = quantizer(input2,q=q,quantizerType='odd',levels=levels)
    evenThreeLevelsInput1 = quantizer(input1,q=q,quantizerType='even',levels=levels)
    evenThreeLevelsInput2 = quantizer(input2,q=q,quantizerType='even',levels=levels)
    
    oddNoise1 = quantizer(0.1*input1,q=q,quantizerType='even',levels=levels)
    oddNoise2 = quantizer(input2,q=q,quantizerType='even',levels=levels)
    oddNNcovs[j] = NP.mean(oddNoise1*oddNoise2)

    autocorrs[j] = NP.mean(input1*input1)
    oddThreeAutocorrs[j] = NP.mean(oddThreeLevelsInput1*oddThreeLevelsInput1)
    evenThreeAutocorrs[j] = NP.mean(evenThreeLevelsInput1*evenThreeLevelsInput1)

    covs[j] = NP.mean(input1*input2)
    threeCovs[j] = NP.mean(oddThreeLevelsInput1*oddThreeLevelsInput2)

    rhos[j] = covs[j]/(NP.std(input1)*NP.std(input2))
    threeRhos[j] = threeCovs[j]/(NP.std(oddThreeLevelsInput1)*NP.std(oddThreeLevelsInput2))
        
PL.figure()
PL.plot(alphas,rhos,'.',label='rhos')
PL.plot(alphas,threeRhos,'.',label='threeRhos')
PL.xlabel('input power')
PL.ylabel('output cross-correlation')
PL.grid()
PL.legend()

PL.figure()
PL.plot(alphas,covs,'.',label='covs')
PL.plot(alphas,threeCovs,'.',label='threeCovs')
PL.xlabel('input power')
PL.ylabel('output cross-power')
PL.grid()
PL.legend()

PL.figure()
PL.plot(alphas,(autocorrs - autocorrs.min())/(autocorrs.max() - autocorrs.min()),'.',label='autocorrs')
PL.plot(alphas,(oddThreeAutocorrs - oddThreeAutocorrs.min() + 1e-6)/(oddThreeAutocorrs.max() - oddThreeAutocorrs.min() + 1e-6),'.',label='oddAutocorrs')
PL.plot(alphas,(evenThreeAutocorrs - evenThreeAutocorrs.min() + 1e-6)/(evenThreeAutocorrs.max() - evenThreeAutocorrs.min() + 1e-6),'.',label='evenAutocorrs')
PL.xlabel('input power')
PL.ylabel('output power')
PL.title(('Levels = %d, q = %d, input ' + r'$\sigma$ %.1f')%(levels, q, input1.std()))
PL.ylim(0.,1.5)
PL.grid()
PL.legend()
