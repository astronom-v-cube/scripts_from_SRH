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

    
L = 500
N = 500
M = 3
scaling = 1000
q = 500
levels = 7

alphas = NP.linspace(0,10,L)

covs = NP.zeros((M,L))
oddAutoCorr1 = NP.zeros((M,L))
oddAutoCorr2 = NP.zeros((M,L))
oddCovs = NP.zeros((M,L))
evenAutoCorr1 = NP.zeros((M,L))
evenAutoCorr2 = NP.zeros((M,L))
evenCovs = NP.zeros((M,L))

for i in range(M):
    diffSigma = .7 - i*.2
    for j in range(L):
        alpha = alphas[j]
        inputSignal = (NP.random.randn(N)*scaling).astype(int)
        noise1 = (NP.random.randn(N)*scaling).astype(int)
        noise2 = (NP.random.randn(N)*scaling).astype(int)
        input1 = diffSigma*(NP.sqrt(alpha)*inputSignal + noise1).astype(int)
        input2 = (NP.sqrt(alpha)*inputSignal + noise2).astype(int)
        
        oddInput1 = quantizer(input1,q=q,quantizerType='odd',levels=levels)
        evenInput1 = quantizer(input1,q=q,quantizerType='even',levels=levels)
        oddInput2 = quantizer(input2,q=q,quantizerType='odd',levels=levels)
        evenInput2 = quantizer(input2,q=q,quantizerType='even',levels=levels)

        covs[i,j] = NP.mean(input1*input2)
        oddAutoCorr1[i,j] = NP.mean(oddInput1*oddInput1)
        oddAutoCorr2[i,j] = NP.mean(oddInput2*oddInput2)
        oddCovs[i,j] = NP.mean(oddInput1*oddInput2)
        evenCovs[i,j] = NP.mean(evenInput1*evenInput2)
        evenAutoCorr1[i,j] = NP.mean(evenInput1*evenInput1)
        evenAutoCorr2[i,j] = NP.mean(evenInput2*evenInput2)
    
fig = PL.figure()
fig.suptitle(r'$\sigma = %d, $q = %d'%(scaling,q))
pl = fig.subplots(nrows = M, ncols = 3)
for i in range(M):
    diffSigma = .7 - i*.2
    pl[i,0].plot(alphas,oddAutoCorr1[i],'.',label='O(%.2f I1 I1)'%diffSigma)
    pl[i,0].plot(alphas,evenAutoCorr1[i],'.',label='E(%.2f I1 I1)'%diffSigma)
    pl[i,1].plot(alphas,oddCovs[i],'.',label='O(%.2f I1) O(I2)'%diffSigma)
    pl[i,1].plot(alphas,evenCovs[i],'.',label='E(%.2f I1) E(I2)'%diffSigma)
    pl[i,2].plot(alphas,oddAutoCorr2[i],'.',label='O(I2 I2)')
    pl[i,2].plot(alphas,evenAutoCorr2[i],'.',label='E(I2 I2)')

    if i == M-1:
        pl[i,0].set_xlabel('input power')
        pl[i,1].set_xlabel('input power')
        pl[i,2].set_xlabel('input power')
    pl[i,0].set_ylabel('output power')
    pl[i,1].set_ylabel('output power')
    pl[i,2].set_ylabel('output power')
    pl[i,0].grid()
    pl[i,1].grid()
    pl[i,2].grid()

    pl[i,0].legend()
    pl[i,1].legend()
    pl[i,2].legend()
    
