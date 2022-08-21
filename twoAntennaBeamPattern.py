#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 08:15:50 2020

@author: svlesovoi
"""
import numpy as NP;
import pylab as PL;
import srhArray

class SrhTwoAntennaPattern():
    def onCalc(self, hourAngle, declination, frequency, diameter, antennaA, antennaB):
        self.fillDish(diameter, frequency)
        
        uvw = self.SRH.baseline2uvw(hourAngle, declination, antennaA, antennaB) / (3e8 / frequency)
        u0 = self.M//2
        v0 = self.M//2
        u1 = self.M//2 +int(uvw[1] / self.uvPerPix + .5)
        v1 = self.M//2 -int(uvw[0] / self.uvPerPix + .5)
        
        self.uvPlainSP[:,:] = complex(0,0)
        self.uvPlainCP[:,:] = complex(0,0)
        self.uvPlainSP[u0 - self.N//2:u0 + self.N//2,v0 - self.N//2:v0 + self.N//2] = self.dish*complex(0,-1)
        self.uvPlainSP[u1 - self.N//2:u1 + self.N//2,v1 - self.N//2:v1 + self.N//2] += self.dish
        self.uvPlainCP[u0 - self.N//2:u0 + self.N//2,v0 - self.N//2:v0 + self.N//2] = self.dish
        self.uvPlainCP[u1 - self.N//2:u1 + self.N//2,v1 - self.N//2:v1 + self.N//2] += self.dish*complex(0,-1)
        self.beamPattern = NP.flipud(self.fftBeam(self.fftConvolution(self.uvPlainSP, self.uvPlainCP).real).real)
        
    def fillDish(self, diameter, frequency):
        radius = diameter / 2 / (3e8 / frequency) / self.uvPerPix
        for i in range(self.N):
            for j in range(self.N):
                x=i - self.N/2
                y=j - self.N/2
                self.dish[i, j] = 0.
                if (NP.sqrt(x**2 + y**2) < radius):
                    self.dish[i, j] = 1.
                    if (x > 0):
                        self.shadowedDish[i, j] = 1.

    def shift2D(self, arr):
        return NP.roll(arr, (arr.shape[0]//2, arr.shape[0]//2), axis=(0,1))
                       
    def fftConvolution(self, arr1, arr2):
        size = arr1.shape[0]
        return NP.roll((NP.fft.ifft2(NP.fft.fft2(arr1) * NP.conjugate(NP.fft.fft2(arr2)))),(size//2,size//2),axis=(0,1))
        
    def fftBeam(self, uvArr):
        return self.shift2D(NP.fft.fft2(self.shift2D(uvArr)))
        
    def __init__(self, fovPixels, arcSecPerPixel):
        self.N = 128
        self.M = fovPixels
        self.meterPerPix = 0.05
        self.pixPerMeter = int(1 / self.meterPerPix + .5)
        self.uvPerPix = 1.
        self.dish = NP.zeros((self.N, self.N))
        self.shadowedDish = NP.zeros((self.N, self.N))
        self.SRH = srhArray.SrhArray()
        self.uvPlainSP = NP.zeros((self.M, self.M),dtype='complex')
        self.uvPlainCP = NP.zeros((self.M, self.M),dtype='complex')
        self.beamPattern = NP.zeros((self.M, self.M))

srhTwoAntennaBeam = SrhTwoAntennaPattern(1024,4.9)
srhTwoAntennaBeam.onCalc(0.1, 0, 5e3*1e6, 1.8, 64, 65)
PL.imshow(srhTwoAntennaBeam.beamPattern)


