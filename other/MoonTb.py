#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  1 04:20:33 2021

@author: svlesovoi

https://www.govinfo.gov/content/pkg/GOVPUB-C13-16b8b2d0e4f8fd358435e9c9ce0353e6/pdf/GOVPUB-C13-16b8b2d0e4f8fd358435e9c9ce0353e6.pdf

"""

import numpy as NP
from scipy import constants
import scipy.optimize as opt

def fitFunc(f, A, B, C, D):
#    return A + B*f + C*f**-3.
    return A + B*f + C*f**D

class MoonTb():
    def __init__(self):
        self.frequency = NP.array(  [0.6,   0.86,   1.5,    3.13,   9.375,  18.75,  37.5,   75.0])# frequency [GHz]
        self.Tb = NP.array(         [.240,  .235,   .225,   .218,   .210,   .207,   .205,   .203])# brightness temperature [1e3K]
        self.guess = [.21,0.01,0.01,-3]
        self.fitTbParams, _ = opt.curve_fit(fitFunc, self.frequency, self.Tb, p0=self.guess)
        self.moonDiskRadius = NP.deg2rad(1000/3600)
        
    def getTbAtFrequency(self, f):
        return fitFunc(f, self.fitTbParams[0],self.fitTbParams[1],self.fitTbParams[2],self.fitTbParams[3])

    def getSfuAtFrequency(self, f):
        return 2*constants.k*self.getTbAtFrequency(f)*1e3/(constants.c/(f*1e9))**2 * NP.pi*self.moonDiskRadius**2 / 1e-22

        