#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 00:20:49 2020

@author: svlesovoi
"""
import os, fnmatch;
import numpy as NP;
import pylab as PL;
from astropy.io import fits
from BadaryRAO import BadaryRAO
from sunpy import coordinates
from matplotlib.ticker import (MultipleLocator)
import skimage
from skimage import filters
from matplotlib.patches import Ellipse

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename))
    result.sort()
    return result

useFittedData = True
showScattered = True
showEllipses = False
showEUV = False

diameterPath = '/home/svlesovoi/Documents/Python Scripts/srhDiameterAngle'
#diameterPath = '/home/sergeylesovoi/Python Scripts/srhDiameterAngle'

fitNames = findFits(diameterPath, '*2020*.fits')
ephDate = '2020-05-18'
RAO = BadaryRAO(ephDate)
solarOptDiameter = NP.rad2deg(RAO.sunObject.radius)*60*2
solarDiameter = 13.927e8
sunRadiusConst = 16. #arcmin

diameterDates = []
diameterFreq = []
sunPAngle = []
sunRadius = []
leftDiameterVsFreq = []
rightDiameterVsFreq = []
leftHeightAbovePhtosphere = []
rightHeightAbovePhtosphere = []
leftDiameterMean = []
rightDiameterMean = []
leftHeightMean = []
rightHeightMean = []
leftBeamAngle = []
rightBeamAngle = []

dRScaledF = 1.022211077227149e-11 #[arcsec/arcsec / Hz]
RScale0 = 0.95
#radiusScale = NP.linspace(0,1,32)*0.045 + 0.95
#radiusScale = NP.linspace(0,31,32)*(4e9/32)*dRScaledF + RScale0
radiusScale = NP.linspace(0.940,.990,32)

for fitName in fitNames[0:]:
    hdu = fits.open(fitName)
    ephDate = hdu[0].header['DATE-OBS']
    diameterDates.append(ephDate)
    RAO = BadaryRAO(ephDate)
    sunPAngle.append(coordinates.sun.P(ephDate).to_value())
    solarOptDiameter = NP.rad2deg(RAO.sunObject.radius)*60*2
    sunRadius.append(solarOptDiameter / 2)
    diameterFreq.append(hdu[1].data['frequencies'])
    leftBeamAngle.append(hdu[1].data['beam_angle_left'])
    rightBeamAngle.append(hdu[1].data['beam_angle_right'])
    leftDiam = radiusScale*hdu[1].data['diameter_left']
    rightDiam = radiusScale*hdu[1].data['diameter_right']
    
    leftDiameterVsFreq.append(leftDiam)
    rightDiameterVsFreq.append(rightDiam)
    leftDiameterMean.append(leftDiam.mean())
    rightDiameterMean.append(rightDiam.mean())
    leftHeightAbovePhtosphere.append(0.5*((leftDiam / solarOptDiameter) * solarDiameter - solarDiameter))
    rightHeightAbovePhtosphere.append(0.5*((rightDiam / solarOptDiameter) * solarDiameter - solarDiameter))
    leftHeightMean.append((leftDiam / solarOptDiameter).mean())
    rightHeightMean.append((rightDiam / solarOptDiameter).mean())
    leftDiam = 0
    rightDiam = 0

sunPAngle = NP.deg2rad(sunPAngle)

leftHeightAbovePhtosphere = NP.array(leftHeightAbovePhtosphere)
rightHeightAbovePhtosphere = NP.array(rightHeightAbovePhtosphere)

leftHeightAbovePhtosphereMean = NP.mean(leftHeightAbovePhtosphere,0)
rightHeightAbovePhtosphereMean = NP.mean(rightHeightAbovePhtosphere,0)
leftHeightAbovePhtosphereStd = NP.std(leftHeightAbovePhtosphere,0)
rightHeightAbovePhtosphereStd = NP.std(rightHeightAbovePhtosphere,0)
leftHeightMean = NP.mean(leftHeightAbovePhtosphere,1)
rightHeightMean = NP.mean(rightHeightAbovePhtosphere,1)

leftFreqs = NP.array(diameterFreq).reshape(len(diameterFreq)*diameterFreq[0].shape[0])
rightFreqs = NP.array(diameterFreq).reshape(len(diameterFreq)*diameterFreq[0].shape[0])
leftHeights = leftHeightAbovePhtosphere.reshape(len(diameterFreq)*diameterFreq[0].shape[0])
rightHeights = rightHeightAbovePhtosphere.reshape(len(diameterFreq)*diameterFreq[0].shape[0])
commonFreqs = NP.concatenate((leftFreqs,rightFreqs))
commonHeights = NP.concatenate((leftHeights, rightHeights))
commonHeightFunc = NP.poly1d(NP.polyfit(commonFreqs, commonHeights, 2))
leftHeightFunc = NP.poly1d(NP.polyfit(leftFreqs,leftHeights,2))
rightHeightFunc = NP.poly1d(NP.polyfit(rightFreqs,rightHeights,2))
freqs = NP.linspace(3e3,9e3,100)

MenezesFreqs = [3000, 5000, 9000, 11000, 13000, 16000, 17000, 22000, 25000]
MenezesHeights = [8e7, 4.4e7, 2.2e7, 2.3e7, 3.0e7, 2.2e7, 1.23e7, 1.6e7, 1.4e7]
MenezesFunc = NP.poly1d(NP.polyfit(MenezesFreqs, MenezesHeights,4))
MenezesFreqs = NP.linspace(3e3,25e3,100)
    
grayColors = PL.get_cmap('gray')

f193 = fits.open('/home/svlesovoi/sunpy/data/aia' + '20200710' + '_000300_0193.fits')
#f193 = fits.open('/home/svlesovoi/sunpy/data/aia_lev1_1700a_2020_07_10t00_00_52_72z_image_lev1.fits')
f193_x_size = f193[0].header['NAXIS1'] * f193[0].header['CDELT1'] / 60
f193_y_size = f193[0].header['NAXIS2'] * f193[0].header['CDELT2'] / 60

PL.figure()
#if showScattered:
#    for i in range(len(diameterFreq)):
#        PL.plot(diameterFreq[i], leftHeightAbovePhtosphere[i], '.', color=grayColors(255 - int(i/len(diameterFreq)*255)))
#        PL.plot(diameterFreq[i], rightHeightAbovePhtosphere[i], '.', color=grayColors(255- int(i/len(diameterFreq)*255)))
#PL.plot(freqs, leftHeightFunc(freqs))
#PL.plot(freqs, rightHeightFunc(freqs))
#PL.plot(freqs, commonHeightFunc(freqs))

#PL.plot(diameterFreq[0], leftHeightAbovePhtosphereMean,'D')
#PL.plot(diameterFreq[0], rightHeightAbovePhtosphereMean,'D')
#PL.plot(diameterFreq[-1], (leftHeightAbovePhtosphereMean + rightHeightAbovePhtosphereMean)*.5,'D')
PL.errorbar(diameterFreq[-1], (leftHeightAbovePhtosphereMean + rightHeightAbovePhtosphereMean)*.5, fmt='D', yerr = (leftHeightAbovePhtosphereStd + rightHeightAbovePhtosphereStd)*.5)
PL.plot(MenezesFreqs, MenezesFunc(MenezesFreqs))
PL.errorbar(3000, 8e7, fmt='o', yerr=1.2e7)    
PL.errorbar(5000, 4.4e7, fmt='o', yerr=6e6)
PL.errorbar(9000, 2.2e7, fmt='o', yerr=1e6)
PL.errorbar(11000, 2.3e7, fmt='o', yerr=4e6)
PL.errorbar(13000, 3.0e7, fmt='o', yerr=1e6)
PL.errorbar(16000, 2.2e7, fmt='o', yerr=3e6)
PL.errorbar(17000, 1.23e7, fmt='o', yerr=1.1e6)
PL.errorbar(22000, 1.6e7, fmt='o', yerr=.6e6)
PL.errorbar(25000, 1.4e7, fmt='o', yerr=3e6)

PL.grid()
PL.ylim(0,1e8)
PL.xlabel('MHz')
PL.ylabel('meters')
PL.title('Height above phtotsphere')

fig, pl = PL.subplots()
pl.set_xlim(-20,20)
pl.set_ylim(-20,20)
pl.xaxis.set_major_locator(MultipleLocator(5.))
pl.yaxis.set_major_locator(MultipleLocator(5.))
pl.grid()
pl.set_xlabel('arcmin')
pl.set_ylabel('arcmin')
theta = NP.deg2rad(NP.linspace(0,360,360))
radiusColors = PL.get_cmap('rainbow')
lowFreqRadius = []
highFreqRadius = []
lowFreqAngle = []
highFreqAngle = []
leftBeamAngle = NP.array(leftBeamAngle)
rightBeamAngle = NP.array(rightBeamAngle)
lRadius = []
rRadius = []
lAngle = []
rAngle = []

if showEUV:
    pl.imshow(f193[1].data*.1,vmax=205,extent=[-f193_x_size/2,f193_x_size/2,f193_y_size/2,-f193_y_size/2])

for i in range(len(diameterFreq)):
    freqColors = radiusColors(NP.linspace(1.,0.,len(diameterFreq[i])))

    leftR = leftDiameterVsFreq[i]/2. * sunRadiusConst / sunRadius[i]
    if useFittedData:
        leftRadiusFunc = NP.poly1d(NP.polyfit(diameterFreq[i],leftR,2))
        leftR = leftRadiusFunc(diameterFreq[i])

    lRadius.append(leftR)
    lowFreqRadius.append(leftR[0])
    highFreqRadius.append(leftR[-1])
    leftAngle = leftBeamAngle[i] + sunPAngle[i]
    lAngle.append(leftAngle)
    lowFreqAngle.append(leftAngle[0])
    highFreqAngle.append(leftAngle[-1])
    if showScattered:
        pl.scatter(-leftR*NP.sin(leftAngle), -leftR*NP.cos(leftAngle), marker='.', color=freqColors, s=5.5)
        pl.scatter( leftR*NP.sin(leftAngle),  leftR*NP.cos(leftAngle), marker='.', color=freqColors, s=5.5)

    rightR = rightDiameterVsFreq[i]/2. * sunRadiusConst / sunRadius[i]
    if useFittedData:
        rightRadiusFunc = NP.poly1d(NP.polyfit(diameterFreq[i],rightR,2))
        rightR = rightRadiusFunc(diameterFreq[i])
        
    rRadius.append(rightR)
    lowFreqRadius.append(rightR[0])
    highFreqRadius.append(rightR[-1])
    rightAngle = rightBeamAngle[i] - sunPAngle[i]
    rAngle.append(rightAngle)
    lowFreqAngle.append(rightAngle[0])
    highFreqAngle.append(rightAngle[-1])
    if showScattered:
        pl.scatter(-rightR*NP.sin(rightAngle), -rightR*NP.cos(rightAngle), marker='.', color=freqColors, s=1.5)
        pl.scatter( rightR*NP.sin(rightAngle),  rightR*NP.cos(rightAngle), marker='.', color=freqColors, s=1.5)

lowFreqRadiusMean = NP.mean(lowFreqRadius)
highFreqRadiusMean = NP.mean(highFreqRadius)
pl.plot(sunRadiusConst*NP.cos(theta), sunRadiusConst*NP.sin(theta), color='orange', label='Solar disk')
#----------------------------------------------------
#lowFreqLeftX = -NP.array(lowFreqRadius)*NP.sin(NP.array(lowFreqAngle))
#lowFreqLeftY = -NP.array(lowFreqRadius)*NP.cos(NP.array(lowFreqAngle))
#lowFreqLeftX = NP.concatenate((lowFreqLeftX, NP.array(lowFreqRadius)*NP.sin(NP.array(lowFreqAngle))))
#lowFreqLeftY = NP.concatenate((lowFreqLeftY, NP.array(lowFreqRadius)*NP.cos(NP.array(lowFreqAngle))))
#lowFreqXY = NP.zeros((lowFreqLeftX.shape[0],2))
#lowFreqXY[:,0] = lowFreqLeftX
#lowFreqXY[:,1] = lowFreqLeftY
#
#highFreqLeftX = -NP.array(highFreqRadius)*NP.sin(NP.array(highFreqAngle))
#highFreqLeftY = -NP.array(highFreqRadius)*NP.cos(NP.array(highFreqAngle))
#highFreqLeftX = NP.concatenate((highFreqLeftX, NP.array(highFreqRadius)*NP.sin(NP.array(highFreqAngle))))
#highFreqLeftY = NP.concatenate((highFreqLeftY, NP.array(highFreqRadius)*NP.cos(NP.array(highFreqAngle))))
#highFreqXY = NP.zeros((highFreqLeftX.shape[0],2))
#highFreqXY[:,0] = highFreqLeftX
#highFreqXY[:,1] = highFreqLeftY
#
#sunEll = skimage.measure.EllipseModel()
#
#sunEll.estimate(lowFreqXY)
#xc, yc, a, b, theta = sunEll.params
#ell_patch = Ellipse((xc, yc), 2*a, 2*b, theta*180/NP.pi, edgecolor=freqColors[0], facecolor='none', label='%d MHz, (%.2f, %.2f)'%(diameterFreq[0][0], a, b), linewidth=.3)
#pl.add_patch(ell_patch)
#
#sunEll.estimate(highFreqXY)
#xc, yc, a, b, theta = sunEll.params
#ell_patch = Ellipse((xc, yc), 2*a, 2*b, theta*180/NP.pi, edgecolor=freqColors[-1], facecolor='none',label='%d MHz, (%.2f, %.2f)'%(diameterFreq[0][-1], a, b), linewidth=.3)
#pl.add_patch(ell_patch)
#----------------------------------------------------
if showEllipses:
    lRadius = NP.array(lRadius)
    lAngle = NP.array(lAngle)
    rRadius = NP.array(rRadius)
    rAngle = NP.array(rAngle)
    
    sunEll = skimage.measure.EllipseModel()
    xRadius = []
    yRadius = []
    for i in range(lRadius.shape[1]):
        ellX = -lRadius[:,i]*NP.sin(lAngle[:,i])
        ellX = NP.concatenate((ellX, lRadius[:,i]*NP.sin(lAngle[:,i])))
        ellX = NP.concatenate((ellX,-rRadius[:,i]*NP.sin(rAngle[:,i])))
        ellX = NP.concatenate((ellX, rRadius[:,i]*NP.sin(rAngle[:,i])))
     
        ellY = -lRadius[:,i]*NP.cos(lAngle[:,i])
        ellY = NP.concatenate((ellY, lRadius[:,i]*NP.cos(lAngle[:,i])))
        ellY = NP.concatenate((ellY,-rRadius[:,i]*NP.cos(rAngle[:,i])))
        ellY = NP.concatenate((ellY, rRadius[:,i]*NP.cos(rAngle[:,i])))
        
        ellXY = NP.zeros((ellX.shape[0],2))
        ellXY[:,0] = ellX
        ellXY[:,1] = ellY
        sunEll.estimate(ellXY)
        yc, xc, a, b, theta = sunEll.params
        xRadius.append(a)
        yRadius.append(b)
        if i == 0 or i == lRadius.shape[1] - 1:
            ell_patch = Ellipse((xc, yc), 2*a, 2*b, theta*180/NP.pi, edgecolor=freqColors[i], facecolor='none', label='%d MHz, Rx %.2f, Ry %.2f'%(diameterFreq[-1][i], a, b), linewidth=.3)
        else:
            ell_patch = Ellipse((xc, yc), 2*a, 2*b, theta*180/NP.pi, edgecolor=freqColors[i], facecolor='none', linewidth=.3)
        pl.add_patch(ell_patch)
#---------------------------------------------------
pl.legend()
if useFittedData:
    kindData = 'fitted'
else:
    kindData = 'data'
pl.set_title('Solar radii %s, %s - %s'%(kindData, diameterDates[0], diameterDates[-1]))
