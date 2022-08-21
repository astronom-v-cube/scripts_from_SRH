#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 08:15:50 2020

@author: svlesovoi
"""
import numpy as NP;
import pylab as PL;
from astropy.coordinates import earth
from astropy.io import fits
from BadaryRAO import BadaryRAO
import scipy.special
from sunpy import coordinates
import srhArray
from scipy import ndimage

srhData = True
EUVdata = True
pAngleRotate = True
useFittedVisibility = True
maxDiameter = 39
minDiameter = 26
ephDate = '2020-09-13'
vis0 = 60
vis1 = 90
dt_major = 3600.
dt_minor = 900.
J_pi = scipy.special.jn_zeros(1,1)[0] / NP.pi
leftModelDia = []
rightModelDia = []
leftHA = []
rightHA = []
leftGp = []
rightGp = []
minFreq = []

class qSunModel():
    def __init__(self, N, radius, value):
        self.disk = NP.zeros((N,N))
        for i in range(N):
            for j in range(N):
                x = i - N/2
                y = j - N/2
                if (NP.sqrt((1.1*x)**2 + y**2) < radius):
                    self.disk[i, j] = value
        self.disk = scipy.ndimage.rotate(self.disk, -25, reshape=False)
    
class SrhTwoAntennaPattern():
    def synthPattern(self, hourAngle, declination, frequency, diameter, antennaA, antennaB):
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
        
    def calcPattern(self, hourAngle, declination, frequency, diameter, antennaA, antennaB, FOV, size):
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

earth_W = NP.rad2deg(earth.OMEGA_EARTH.to_value())*3600

fitsDate = ''.join(ephDate.split('-'))
redVis = fits.open('srh_redVis_' + fitsDate + '.fits')
#f193 = fits.open('/home/sergeylesovoi/sunpy/data/aia' + fitsDate + '_000300_0193.fits')
#f193 = fits.open('/home/svlesovoi/sunpy/data/aia' + fitsDate + '_000300_0193.fits')
f193 = fits.open('/home/svlesovoi/sunpy/data/aia' + '20200710' + '_000300_0193.fits')
#f193 = fits.open('/home/sergeylesovoi/sunpy/data/aia' + '20200720_000300_0193.fits')
x_size = f193[0].header['NAXIS1'] * f193[0].header['CDELT1']
y_size = f193[0].header['NAXIS2'] * f193[0].header['CDELT2']

RAO = BadaryRAO(ephDate)
pAngle = NP.deg2rad(coordinates.sun.P(ephDate).to_value())
qSun = qSunModel(f193[0].header['NAXIS1'], int(17.*60/f193[0].header['CDELT1']), 100.)

redFreq = redVis[1].data['frequencies']
redTime = redVis[2].data['time']
redA = redVis[3].data['antennaA']
redB = redVis[3].data['antennaB']
redRcp = redVis[4].data['VIS_RCP']
redLcp = redVis[4].data['VIS_LCP']
redRcpArr = redRcp.reshape(redFreq.shape[0],redTime.shape[1],redA.shape[0])
redLcpArr = redLcp.reshape(redFreq.shape[0],redTime.shape[1],redA.shape[0])
red2RcpAbs = NP.mean(NP.abs(redRcpArr[:,:,vis0:vis1]),2)
red2LcpAbs = NP.mean(NP.abs(redLcpArr[:,:,vis0:vis1]),2)
red2RcpPha = NP.angle(redRcpArr[:,:,vis0:vis1])
red2LcpPha = NP.angle(redLcpArr[:,:,vis0:vis1])
for f in range(redFreq.shape[0]):
    for v in range(vis1 - vis0):
        red2RcpPha[f,:,v] = NP.unwrap(red2RcpPha[f,:,v])
        red2LcpPha[f,:,v] = NP.unwrap(red2LcpPha[f,:,v])
red2RcpPha = NP.mean(red2RcpPha,2)
red2LcpPha = NP.mean(red2LcpPha,2)
tempLeftTheta = []
tempLeftVisMinimum = []
tempLeftVisFitted = []
tempRightTheta = []
tempRightVisMinimum = []
tempRightVisFitted = []

if EUVdata:
    if pAngleRotate:
        f193data = ndimage.rotate(NP.flipud(f193[0].data),NP.rad2deg(pAngle),reshape=False, mode='nearest')
    else:
        f193data = f193[0].data
else:
    f193data = qSun.disk
#    qDiskSrh = fits.open('20200510_aver_equatorial.fits')
#    f193data = NP.flipud(scipy.ndimage.zoom(qDiskSrh[0].data,4))
    
    
srhTwoAntennaBeam = SrhTwoAntennaPattern(f193[0].header['NAXIS1'],f193[0].header['CDELT1'])

for f in NP.linspace(0,31,32,dtype='int'):
    hA = NP.deg2rad(earth_W*(redTime[f,:] - RAO.culmination)/3600)
    cosP = NP.sin(hA) * NP.cos(RAO.declination)
    theta = NP.rad2deg(3e8/(redFreq[f]*1e6)/(9.8*NP.sqrt(1 - cosP*cosP)))*60*NP.sign(hA)#lambda / (b sin(p)) signum(h) [arcmin]
    gP =  NP.arctan(NP.tan(hA)*NP.sin(RAO.declination)) + NP.pi/2
    gQ =  NP.arctan(-(NP.sin(RAO.declination)/NP.tan(hA) + NP.cos(RAO.declination)/(NP.sin(hA)*NP.tan(RAO.observatory.lat)))) + NP.pi/2
    gQ[NP.where(hA > 0.)] = NP.pi + gQ[NP.where(hA > 0.)]
    minFreq.append(redFreq[f])

    ind = NP.where((theta > -maxDiameter/J_pi) * (theta < -minDiameter/J_pi))[0]
    if srhData:
        modelVisibility = red2RcpAbs[f,ind] + red2LcpAbs[f,ind]
        tempLeftVisMinimum.append(modelVisibility)
    else:
        modelVisibility = []
        for t in ind:
            srhTwoAntennaBeam.calcPattern(hA[t], RAO.declination, redFreq[f]*1e6, 1.8, 64, 66, NP.deg2rad(x_size/3600), f193[0].header['NAXIS1'])
            modelVisibility.append(NP.abs(NP.mean(srhTwoAntennaBeam.beamPattern*f193data)))
    if useFittedVisibility:
        visibilityFunc = NP.poly1d(NP.polyfit(theta[ind],modelVisibility,4))
        fittedVisibility = visibilityFunc(theta[ind])
        tempLeftVisFitted.append(fittedVisibility)
        tempLeftTheta.append(theta[ind])
        leftModelDia.append(-J_pi*theta[ind[NP.argmin(fittedVisibility)]])#1.22 lambda / ((b sin(p)) signum(h))
        leftGp.append(gP[ind[NP.argmin(fittedVisibility)]])
    else:
        leftModelDia.append(-J_pi*theta[ind[NP.argmin(modelVisibility)]])
        leftGp.append(gP[ind[NP.argmin(modelVisibility)]])

    ind = NP.where((theta < maxDiameter/J_pi) * (theta > minDiameter/J_pi))[0]
    if srhData:
        modelVisibility = red2RcpAbs[f,ind] + red2LcpAbs[f,ind]
        tempRightVisMinimum.append(modelVisibility)
    else:
        modelVisibility = []
        for t in ind:
            srhTwoAntennaBeam.calcPattern(hA[t], RAO.declination, redFreq[f]*1e6, 1.8, 64, 66, NP.deg2rad(x_size/3600), f193[0].header['NAXIS1'])
            modelVisibility.append(NP.abs(NP.mean(srhTwoAntennaBeam.beamPattern*f193data)))
    if useFittedVisibility:
        visibilityFunc = NP.poly1d(NP.polyfit(theta[ind],modelVisibility,4))
        fittedVisibility = visibilityFunc(theta[ind])
        tempRightTheta.append(theta[ind])
        tempRightVisFitted.append(fittedVisibility)
        rightModelDia.append(J_pi*theta[ind[NP.argmin(fittedVisibility)]])
        rightGp.append(gP[ind[NP.argmin(fittedVisibility)]])
    else:
        rightModelDia.append(J_pi*theta[ind[NP.argmin(modelVisibility)]])
        rightGp.append(gP[ind[NP.argmin(modelVisibility)]])

leftDiaFunc = NP.poly1d(NP.polyfit(minFreq,leftModelDia,2))
rightDiaFunc = NP.poly1d(NP.polyfit(minFreq,rightModelDia,2))

PL.figure(figsize = (6,3))
PL.title(ephDate)
PL.plot(minFreq,leftModelDia, '.')
PL.plot(minFreq,rightModelDia, '.')
PL.plot(minFreq,leftDiaFunc(minFreq))
PL.plot(minFreq,rightDiaFunc(minFreq))
PL.ylim(30,40)
PL.ylabel('solar diameter [arcmin]')
PL.xlabel('frequency [MHz]')
PL.grid()

PL.figure()
PL.title(ephDate)
for f in range(len(minFreq)):
    xx =  960*NP.sin(leftGp[f])
    yy = -960*NP.cos(leftGp[f])
    PL.plot([-xx,xx],[-yy,yy],color='green',linewidth=0.3)
    xx =  960*NP.sin(rightGp[f])
    yy = -960*NP.cos(rightGp[f])
    PL.plot([-xx,xx],[-yy,yy],color='red',linewidth=0.3)

alpha = NP.linspace(0,2*NP.pi,360)
PL.plot(16*60 * NP.cos(2*NP.pi*alpha), 16*60 * NP.sin(2*NP.pi*alpha), '.')
PL.plot(17*60 * NP.cos(2*NP.pi*alpha), 17*60 * NP.sin(2*NP.pi*alpha), '.')
PL.plot(18*60 * NP.cos(2*NP.pi*alpha), 18*60 * NP.sin(2*NP.pi*alpha), '.')

uvPerPixel = 1./(NP.deg2rad(x_size * 3600))
PL.imshow(f193data,vmax=100,extent=[-x_size/2,x_size/2,y_size/2,-y_size/2])

if srhData:
    freqsColumn = fits.Column(name='frequencies',format='E',array=NP.array(minFreq))
    leftDiameter = fits.Column(name='diameter_left',format='E',array=NP.array(leftModelDia))
    rightDiameter = fits.Column(name='diameter_right',format='E',array=NP.array(rightModelDia))
    leftBeamAngle = fits.Column(name='beam_angle_left',format='E',array=NP.array(leftGp))
    rightBeamAngle = fits.Column(name='beam_angle_right',format='E',array=NP.array(rightGp))
    
    diameterTableHdu = fits.BinTableHDU.from_columns([freqsColumn, leftDiameter, rightDiameter, leftBeamAngle, rightBeamAngle])
    
    pHeader = fits.Header();
    pHeader['DATE-OBS']     = redVis[0].header['DATE-OBS'];
    pHeader['TIME-OBS']     = redVis[0].header['TIME-OBS'];
    pHeader['INSTRUME']     = redVis[0].header['INSTRUME'];
    pHeader['ORIGIN']       = redVis[0].header['ORIGIN'];
    pHeader['OBS-LAT']      = redVis[0].header['OBS-LAT'];
    pHeader['OBS-LONG']     = redVis[0].header['OBS-LONG'];
    pHeader['OBS-ALT']      = redVis[0].header['OBS-ALT'];
    pHdu = fits.PrimaryHDU(header=pHeader);
    hduList = fits.HDUList([pHdu, diameterTableHdu]);
    hduList.writeto('srh_diameter_' + ephDate + '.fits',clobber=True);
    hduList.close();

del redRcp
del redLcp
del redRcpArr
del redLcpArr
del red2RcpAbs
del red2LcpAbs
del red2RcpPha
del red2LcpPha
del f193
del redVis
