#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 08:23:06 2022

@author: sergeyvlesovoi
"""


import pylab as PL
import numpy as NP
import scipy as SP

def quantizer(inputSignal, q = 1, levels = 1,quantizerType = 'even'):
    if quantizerType == 'even':
        return NP.clip(inputSignal//q*q + q/2,-(levels - .5)*q,(levels - .5)*q)
    else:
        return NP.clip((inputSignal + q/2)//q*q,-levels*q,levels*q)

def quSigma(sigma, nBits):
    maxLevel = 2 ** (nBits-1) - 1
    levels = 0
    for k in range(maxLevel):
        levels += (2*k + 1) * SP.special.erf((k+.5)/(NP.sqrt(2)*sigma))
    return NP.sqrt(maxLevel**2 - levels)

def quSigma3levels(sigma, maxLevel):
    return NP.sqrt(maxLevel**2 * (1 - SP.special.erf(.5/(NP.sqrt(2)*sigma))))

N = 500
L = 500

alphas = NP.linspace(0,10,L)
rhos = NP.zeros(L)
threeRhos = NP.zeros(L)
covs = NP.zeros(L)
threeCovs = NP.zeros(L)
autocorrs = NP.zeros(L)
oddAutocorrs = NP.zeros((5,L))
evenAutocorrs = NP.zeros((5,L))

scaling = 100
levels = 7

fig = PL.figure()
pl = fig.subplots(nrows=1,ncols=5)

for q in range(5):
    qS = scaling*q/2 + .00001
    sigma = 0
    for j in range(L):
        alpha = alphas[j]
        inputSignal = (NP.random.randn(N)*scaling).astype(int)
        noise1 = (NP.random.randn(N)*scaling).astype(int)
        input1 = (NP.sqrt(alpha)*inputSignal + noise1).astype(int)
        sigma = input1.std()

        oddInput1 = quantizer(input1,q=qS,quantizerType='odd',levels=levels)
        evenInput1 = quantizer(input1,q=qS,quantizerType='even',levels=levels)
    
        autocorrs[j] = NP.mean(input1*input1)
        oddAutocorrs[q,j] = NP.mean(oddInput1*oddInput1)
        evenAutocorrs[q,j] = NP.mean(evenInput1*evenInput1)
    
    pl[q].title.set_text(('q = %d, ' + r'$\sigma$ = %.1f')%(qS, sigma))
    pl[q].plot(alphas,(autocorrs - autocorrs.min())/(autocorrs.max() - autocorrs.min()),'.',label='auto')
    pl[q].plot(alphas,(oddAutocorrs[q] - oddAutocorrs[q].min() + 1e-6)/(oddAutocorrs[q].max() - oddAutocorrs[q].min() + 1e-6),'.',label='odd')
    pl[q].plot(alphas,(evenAutocorrs[q] - evenAutocorrs[q].min() + 1e-6)/(evenAutocorrs[q].max() - evenAutocorrs[q].min() + 1e-6),'.',label='even')
#        pl[qSigma].xlabel('input power')
#        pl[qSigma].ylabel('output power')
#        pl[qSigma].title(('Levels = %d, q = %d, input ' + r'$\sigma$ %.1f')%(levels, q, input1.std()))
#        pl[qSigma].ylim(0.,1.5)
    pl[q].grid()
    pl[q].legend()


#sigs = NP.linspace(0,3,1000)
#PL.figure()
#PL.plot(sigs,quSigma(sigs,12),label='2047 levels (12-bit)')
#PL.plot(sigs,quSigma(sigs,4),label='15 levels (4-bit)')
#PL.plot(sigs,quSigma(sigs,2),label='3 levels (2-bit)')
#PL.xlabel(r'$\sigma$')
#PL.ylabel(r'$\hat{\sigma}$')
#PL.legend()
#PL.grid()
#
#PL.figure()
#PL.plot(sigs**2,quSigma(sigs,12)**2,label='2047 levels (12-bit)')
#PL.plot(sigs**2,quSigma(sigs,4)**2,label='15 levels (4-bit)')
#PL.plot(sigs**2,quSigma(sigs,2),label='3 levels (2-bit)')
#PL.xlabel(r'$\sigma^2$')
#PL.ylabel(r'$\hat{\sigma}^2$')
#PL.legend()
#PL.grid()
