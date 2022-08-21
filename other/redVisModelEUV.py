#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 08:15:50 2020

@author: svlesovoi
"""
import numpy as NP;
import pylab as PL;
from astropy.io import fits
from BadaryRAO import BadaryRAO
from scipy.stats import linregress
import scipy.special
from PIL import Image
from sunpy import coordinates
import srhArray
from scipy import ndimage

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
        
    def calcBeam(self, hourAngle, declination, frequency, diameter, antennaA, antennaB, FOV, size):
        uvw = self.SRH.baseline2uvw(hourAngle, declination, antennaA, antennaB) / (3e8 / frequency)
        
        uu = 2*NP.pi*uvw[1]*NP.linspace(-FOV/2,FOV/2,size)
        vv = 2*NP.pi*uvw[0]*NP.linspace(-FOV/2,FOV/2,size)
        ux, vy = NP.meshgrid(uu, vv)
        self.beamPattern = NP.cos(ux + vy)
#        xx = NP.pi*diameter/(3e8/frequency)*NP.linspace(-FOV/2,FOV/2,size)
#        yy = NP.pi*diameter/(3e8/frequency)*NP.linspace(-FOV/2,FOV/2,size)
#        x, y = NP.meshgrid(xx, yy)
#        self.beamPattern *= scipy.special.jv(1, x + y) / (x + y + 0.0001)

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
        self.uvPerPix = 1.#1./NP.deg2rad(fovPixels * arcSecPerPixel / 3600)
        self.dish = NP.zeros((self.N, self.N))
        self.shadowedDish = NP.zeros((self.N, self.N))
        self.SRH = srhArray.SrhArray()
        self.uvPlainSP = NP.zeros((self.M, self.M),dtype='complex')
        self.uvPlainCP = NP.zeros((self.M, self.M),dtype='complex')
        self.beamPattern = NP.zeros((self.M, self.M))

dt_major = 3600.;
dt_minor = 900.;

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

#ephDate = '2020-03-25'
#ephDate = '2020-04-23'
#ephDate = '2020-04-24'
#ephDate = '2020-05-02'
#ephDate = '2020-05-05'
#ephDate = '2020-05-09'
#ephDate = '2020-05-10'
ephDate = '2020-05-11'
RAO = BadaryRAO(ephDate)
pAngle = NP.deg2rad(coordinates.sun.P(ephDate).to_value())

N = 1024
qSun = NP.zeros((N,N))

radius = int(17.*60/2.4)
for i in range(N):
    for j in range(N):
        x = i - N/2
        y = j - N/2
        if (NP.sqrt(x**2 + y**2) < radius):
            qSun[i, j] = 100.

#redVis = fits.open('srh_redVis_20200325.fits')
#redVis = fits.open('srh_redVis_20200423.fits')
#redVis = fits.open('srh_redVis_20200424.fits')
#redVis = fits.open('srh_redVis_20200502.fits')
#redVis = fits.open('srh_redVis_20200505.fits')
#redVis = fits.open('srh_redVis_20200509.fits')
#redVis = fits.open('srh_redVis_20200510.fits')
redVis = fits.open('srh_redVis_20200511.fits')

redFreq = redVis[1].data['frequencies']
redTime = redVis[2].data['time']
redA = redVis[3].data['antennaA']
redB = redVis[3].data['antennaB']
redRcp = redVis[4].data['VIS_RCP']
redLcp = redVis[4].data['VIS_LCP']
redRcpArr = redRcp.reshape(redFreq.shape[0],redTime.shape[1],redA.shape[0])
redLcpArr = redLcp.reshape(redFreq.shape[0],redTime.shape[1],redA.shape[0])
vis0 = 60
vis1 = 90
#red2RcpMean = NP.mean(NP.abs(redRcpArr[:,:,vis0:vis1]),2)
#red2LcpMean = NP.mean(NP.abs(redLcpArr[:,:,vis0:vis1]),2)
red2RcpMean = NP.mean(redRcpArr.real[:,:,vis0:vis1],2)
red2LcpMean = NP.mean(redLcpArr.real[:,:,vis0:vis1],2)

#f193 = fits.open('/home/svlesovoi/sunpy/data/aia20200509_000600_0193.fits')
#f193 = fits.open('/home/svlesovoi/sunpy/data/aia20200510_000300_0193.fits')
f193 = fits.open('/home/svlesovoi/sunpy/data/aia20200511_000300_0193.fits')
x_size = f193[0].header['NAXIS1'] * f193[0].header['CDELT1']
y_size = f193[0].header['NAXIS2'] * f193[0].header['CDELT2']
#f193data = ndimage.rotate(f193[0].data,NP.rad2deg(pAngle),reshape=False)
#f193data = f193[0].data
f193data = qSun
srhTwoAntennaBeam = SrhTwoAntennaPattern(f193[0].header['NAXIS1'],f193[0].header['CDELT1'])

vis = 70
freq = 30
leftModelMinRad = []
rightModelMinRad = []
leftHA = []
rightHA = []
leftGp = []
rightGp = []
minFreq = []
hAbins = 400
modelVisibility = NP.zeros((32,hAbins))
theta = NP.zeros((32,hAbins))
gP = NP.zeros(hAbins)

J_pi = scipy.special.jn_zeros(1,1)[0] / NP.pi

for f in NP.linspace(0,31,32,dtype='int'):
    hA = NP.deg2rad(15*(redTime[f,:] - RAO.culmination)/3600)
    cosP = NP.sin(hA) * NP.cos(RAO.declination)
    minFreq.append(redFreq[f])

    hT = NP.linspace(0,hA.shape[0]-1,hAbins,dtype='int')
    for t in NP.linspace(0,hAbins-1,hAbins,dtype='int'):
        srhTwoAntennaBeam.calcBeam(hA[hT[t]], RAO.declination, redFreq[f]*1e6, 1.8, 64, 66, NP.deg2rad(x_size/3600), f193[0].header['NAXIS1'])
        modelVisibility[f,t] = NP.abs(NP.mean(srhTwoAntennaBeam.beamPattern*f193data))
        theta[f,t] = NP.rad2deg(3e8/(redFreq[f]*1e6)/(9.8*NP.sqrt(1 - cosP[hT[t]]*cosP[hT[t]])))*60*NP.sign(hA[hT[t]])
        gP[t] =  NP.arctan(NP.tan(hA[hT[t]])*NP.sin(RAO.declination)) + NP.pi/2

    ind = NP.where((theta[f] > -32) * (theta[f] < -26))
    subMean = modelVisibility[f][ind]
    subMeanFunc = NP.poly1d(NP.polyfit(theta[f][ind],subMean,3))
    fittedSubMean = subMeanFunc(theta[f][ind])
    leftModelMinRad.append(-J_pi*theta[f][ind[0][NP.argmin(fittedSubMean)]])
    leftGp.append(gP[ind[0][NP.argmin(fittedSubMean)]])

    ind = NP.where((theta[f] < 36) * (theta[f] > 26))
    subMean = modelVisibility[f][ind]
    subMeanFunc = NP.poly1d(NP.polyfit(theta[f][ind],subMean,3))
    fittedSubMean = subMeanFunc(theta[f][ind])
    rightModelMinRad.append(-J_pi*theta[f][ind[0][NP.argmin(fittedSubMean)]])
    rightGp.append(gP[ind[0][NP.argmin(fittedSubMean)]])

fig = PL.figure(figsize = (6,3))
for f in NP.linspace(0,31,32,dtype='int'):
    PL.plot(hA[hT],modelVisibility[f])
    
fig = PL.figure(figsize = (6,3))
PL.plot(redFreq,NP.abs(leftModelMinRad))
PL.plot(redFreq,NP.abs(rightModelMinRad))
PL.ylim(30,40)

PL.figure()
for f in range(len(minFreq)):
    xx =  x_size/2*NP.sin(leftGp[f])
    yy = -y_size/2*NP.cos(leftGp[f])
    PL.plot([-xx,xx],[-yy,yy],color='green',linewidth=0.3)
    xx =  x_size/2*NP.sin(rightGp[f])
    yy = -y_size/2*NP.cos(rightGp[f])
    PL.plot([-xx,xx],[-yy,yy],color='red',linewidth=0.3)

uvPerPixel = 1./(NP.deg2rad(x_size * 3600))
PL.imshow(f193data,vmax=100,extent=[-x_size/2,x_size/2,y_size/2,-y_size/2])

del redVis
del redRcp
del redLcp

del redRcpArr
del redLcpArr
#del red2RcpMean
#del red2LcpMean

