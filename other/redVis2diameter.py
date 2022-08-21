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
import scipy.special
from sunpy import coordinates
import os, fnmatch;

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

srhData = True
pAngleRotate = True
useFittedVisibility = True
maxDiameter = 38
minDiameter = 26
ephDate = '2020-06-20'
vis0 = 60
vis1 = 90
dt_major = 3600.
dt_minor = 900.
J_pi = scipy.special.jn_zeros(1,1)[0] / NP.pi

fitNames = findFits('/home/svlesovoi/Documents/Python Scripts/srhRedVis/', '*.fits')
fitNames.sort()

for fitName in fitNames:
    print(fitName)
    leftModelDia = []
    rightModelDia = []
    leftHA = []
    rightHA = []
    leftGp = []
    rightGp = []
    minFreq = []
    
    redVis = fits.open(fitName)
    ephDate = redVis[0].header['DATE-OBS'].replace('/','-')
    RAO = BadaryRAO(ephDate)
    pAngle = NP.deg2rad(coordinates.sun.P(ephDate).to_value())
    
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
    
    for f in NP.linspace(0,31,32,dtype='int'):
        hA = NP.deg2rad(15*(redTime[f,:] - RAO.culmination)/3600)
        cosP = NP.sin(hA) * NP.cos(RAO.declination)
        theta = NP.rad2deg(3e8/(redFreq[f]*1e6)/(9.8*NP.sqrt(1 - cosP*cosP)))*60*NP.sign(hA)
        gP =  NP.arctan(NP.tan(hA)*NP.sin(RAO.declination)) + NP.pi/2
        gQ =  NP.arctan(-(NP.sin(RAO.declination)/NP.tan(hA) + NP.cos(RAO.declination)/(NP.sin(hA)*NP.tan(RAO.observatory.lat)))) + NP.pi/2
        gQ[NP.where(hA > 0.)] = NP.pi + gQ[NP.where(hA > 0.)]
        minFreq.append(redFreq[f])
    
        ind = NP.where((theta > -maxDiameter/J_pi) * (theta < -minDiameter/J_pi))[0]
        modelVisibility = red2RcpAbs[f,ind] + red2LcpAbs[f,ind]
        if useFittedVisibility:
            visibilityFunc = NP.poly1d(NP.polyfit(theta[ind],modelVisibility,4))
            fittedVisibility = visibilityFunc(theta[ind])
            leftModelDia.append(-J_pi*theta[ind[NP.argmin(fittedVisibility)]])
            leftGp.append(gP[ind[NP.argmin(fittedVisibility)]])
        else:
            leftModelDia.append(-J_pi*theta[ind[NP.argmin(modelVisibility)]])
            leftGp.append(gP[ind[NP.argmin(modelVisibility)]])
    
        ind = NP.where((theta < maxDiameter/J_pi) * (theta > minDiameter/J_pi))[0]
        modelVisibility = red2RcpAbs[f,ind] + red2LcpAbs[f,ind]
        if useFittedVisibility:
            visibilityFunc = NP.poly1d(NP.polyfit(theta[ind],modelVisibility,4))
            fittedVisibility = visibilityFunc(theta[ind])
            rightModelDia.append(J_pi*theta[ind[NP.argmin(fittedVisibility)]])
            rightGp.append(gP[ind[NP.argmin(fittedVisibility)]])
        else:
            rightModelDia.append(J_pi*theta[ind[NP.argmin(modelVisibility)]])
            rightGp.append(gP[ind[NP.argmin(modelVisibility)]])
    
    leftDiaFunc = NP.poly1d(NP.polyfit(minFreq,leftModelDia,2))
    rightDiaFunc = NP.poly1d(NP.polyfit(minFreq,rightModelDia,2))
    
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
    del redVis

#PL.figure(figsize = (6,3))
#PL.title(ephDate)
#PL.plot(minFreq,leftModelDia, '.')
#PL.plot(minFreq,rightModelDia, '.')
#PL.plot(minFreq,leftDiaFunc(minFreq))
#PL.plot(minFreq,rightDiaFunc(minFreq))
#PL.ylim(30,40)
#PL.ylabel('solar diameter [arcmin]')
#PL.xlabel('frequency [MHz]')
#PL.grid()
#
#PL.figure()
#PL.title(ephDate)
#for f in range(len(minFreq)):
#    xx =  960*NP.sin(leftGp[f])
#    yy = -960*NP.cos(leftGp[f])
#    PL.plot([-xx,xx],[-yy,yy],color='green',linewidth=0.3)
#    xx =  960*NP.sin(rightGp[f])
#    yy = -960*NP.cos(rightGp[f])
#    PL.plot([-xx,xx],[-yy,yy],color='red',linewidth=0.3)
#
#alpha = NP.linspace(0,2*NP.pi,360)
#PL.plot(16*60 * NP.cos(2*NP.pi*alpha), 16*60 * NP.sin(2*NP.pi*alpha), '.')
#PL.plot(17*60 * NP.cos(2*NP.pi*alpha), 17*60 * NP.sin(2*NP.pi*alpha), '.')
#PL.plot(18*60 * NP.cos(2*NP.pi*alpha), 18*60 * NP.sin(2*NP.pi*alpha), '.')

