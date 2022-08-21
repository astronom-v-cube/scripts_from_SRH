#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 09:20:53 2022

@author: sergeyvlesovoi
"""


import numpy as NP
import pylab as PL

N = 100000
C = 0.00

sig0 = NP.random.randn(N)

sig1 = C*sig0 + NP.random.randn(N)
histSig1 = NP.histogram(sig1,bins=1000)

sig2 = C*sig0 + NP.random.randn(N)
histSig2 = NP.histogram(sig2,bins=1000)

accT = 100
bins = NP.arange(-10,10,.01)

pow1 = sig1*sig1
accPow1 = NP.convolve(pow1,NP.ones(accT),mode='valid')/accT
histPow1 = NP.histogram(pow1,bins=bins)
histAccPow1 = NP.histogram(accPow1,bins=bins)

pow12 = sig1*sig2
accPow12 = NP.convolve(pow12,NP.ones(accT),mode='valid')/accT
histPow12 = NP.histogram(pow12,bins=bins)
histAccPow12 = NP.histogram(accPow12,bins=bins)

fig, pl= PL.subplots(nrows=3, ncols=1,figsize=(8,12))

pl[0].plot(C*sig0,label='C*sig0')
pl[0].set_xlabel('time')
pl[0].set_ylabel('sig0')
pl[0].set_ylim(-10,10)
pl[0].grid()
pl[0].legend()

pl[1].plot(sig1,label='sig1')
pl[1].set_xlabel('time')
pl[1].set_ylabel('sig1')
pl[1].set_ylim(-10,10)
pl[1].grid()
pl[1].legend()

pl[2].plot(sig2,label='sig2')
pl[2].set_xlabel('time')
pl[2].set_ylabel('sig2')
pl[2].set_ylim(-10,10)
pl[2].grid()
pl[2].legend()

PL.show()

fig, pl= PL.subplots(nrows=3, ncols=1,figsize=(8,12))
#fig.suptitle('C=%.1f'%C)

pl[0].plot(pow1,label='radiometer',color='red')
pl[0].plot(accPow1,label='accRadiometer',color='orange')
pl[0].set_xlabel('time')
pl[0].set_ylabel('sig1*sig1')
pl[0].set_ylim(-100,100)
pl[0].grid()
pl[0].legend()

pl[1].plot(pow12,label='interferometer',color='blue')
pl[1].plot(accPow12,label='accInterfeormeter',color='magenta')
pl[1].set_xlabel('time')
pl[1].set_ylabel('sig1*sig2')
pl[1].set_ylim(-100,100)
pl[1].grid()
pl[1].legend()

pl[2].plot(histPow1[1][1:],histPow1[0],label='radiometer',color='red')
pl[2].plot(histAccPow1[1][1:],histAccPow1[0],label='accRadiometer %.2f'%accPow1.std(),color='orange')

pl[2].plot(histPow12[1][1:],histPow12[0],label='interfeormeter',color='blue')
pl[2].plot(histAccPow12[1][1:],histAccPow12[0],label='accInterfeormeter %.2f'%accPow12.std(),color='magenta')
pl[2].set_xlim(-7,7)
pl[2].set_ylim(0,5000)
pl[2].set_xlabel('value')
pl[2].set_ylabel('count')
pl[2].grid()
pl[2].legend()
PL.show()
