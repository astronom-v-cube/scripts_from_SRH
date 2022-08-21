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
ephDate = '2020-05-09'
#ephDate = '2020-05-10'
RAO = BadaryRAO(ephDate)
pAngle = NP.deg2rad(coordinates.sun.P(ephDate).to_value())

#redVis = fits.open('srh_redVis_20200325.fits')
#redVis = fits.open('srh_redVis_20200423.fits')
#redVis = fits.open('srh_redVis_20200424.fits')
#redVis = fits.open('srh_redVis_20200502.fits')
#redVis = fits.open('srh_redVis_20200505.fits')
redVis = fits.open('srh_redVis_20200509.fits')
#redVis = fits.open('srh_redVis_20200510.fits')

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
red2RcpMean = NP.mean(NP.abs(redRcpArr[:,:,vis0:vis1]),2)
red2LcpMean = NP.mean(NP.abs(redLcpArr[:,:,vis0:vis1]),2)
red2U = red2RcpMean.copy()

vis0 = 15
vis1 = 46
red1RcpMean = NP.mean(NP.abs(redRcpArr[:,:,vis0:vis1]),2)
red1LcpMean = NP.mean(NP.abs(redLcpArr[:,:,vis0:vis1]),2)
red1RcpMeanPhase = NP.mean(NP.angle(redRcpArr[:,:,vis0:vis0+1]),2)
red1LcpMeanPhase = NP.mean(NP.angle(redLcpArr[:,:,vis0:vis0+1]),2)
red1U = red1RcpMean.copy()

vis = 70
freq = 30
leftMinRad = []
rightMinRad = []
leftGp = []
rightGp = []
minFreq = []
t_cal0 = 90
#t_cal0 = 360
t_cal1 = 2460
dt_cal = 80
J_pi = scipy.special.jn_zeros(1,1)[0] / NP.pi

for f in NP.linspace(0,31,32,dtype='int'):
    hA = NP.deg2rad(15*(redTime[f,:] - RAO.culmination)/3600)
    gP =  NP.arctan(NP.tan(hA)*NP.sin(RAO.declination)) + NP.pi/2
    gQ =  NP.arctan(-(NP.sin(RAO.declination)/NP.tan(hA) + NP.cos(RAO.declination)/(NP.sin(hA)*NP.tan(RAO.observatory.lat)))) + NP.pi/2
    gQ[NP.where(hA > 0.)] = NP.pi + gQ[NP.where(hA > 0.)]
    
    cosP = NP.sin(hA) * NP.cos(RAO.declination)
    theta = NP.rad2deg(3e8/(redFreq[f]*1e6)/(9.8*NP.sqrt(1 - cosP*cosP)))*60*NP.sign(hA)
    red1U[f,:] = 4.9*NP.sqrt(1 - cosP*cosP)*NP.sign(hA)/(3e8/(redFreq[f]*1e6))
    red2U[f,:] = 2*red1U[f,:]
    minFreq.append(redFreq[f])
#    red2RcpMean[f,t_cal0:t_cal0 + dt_cal] = 0.1
#    red2LcpMean[f,t_cal0:t_cal0 + dt_cal] = 0.1
#
#    red2RcpMean[f,t_cal1:t_cal1 + dt_cal] = 0.1
#    red2LcpMean[f,t_cal1:t_cal1 + dt_cal] = 0.1

    ind = NP.where((theta > -30.5) * (theta < -26.5))
    subMean = red2RcpMean[f,:][ind] + red2LcpMean[f,:][ind]
    subMeanFunc = NP.poly1d(NP.polyfit(theta[ind],subMean,3))
    fittedSubMean = subMeanFunc(theta[ind])
    leftMinRad.append(-J_pi*theta[ind[0][NP.argmin(fittedSubMean)]])
    leftGp.append(gP[ind[0][NP.argmin(fittedSubMean)]])

    ind = NP.where((theta < 36) * (theta > 26))
    subMean = red2RcpMean[f,:][ind] + red2LcpMean[f,:][ind]
    subMeanFunc = NP.poly1d(NP.polyfit(theta[ind],subMean,3))
    fittedSubMean = subMeanFunc(theta[ind])
    rightMinRad.append(J_pi*theta[ind[0][NP.argmin(fittedSubMean)]])
    rightGp.append(gP[ind[0][NP.argmin(fittedSubMean)]])

leftRadSlope, leftRadIntercept, rV, pV, leftErr = linregress(minFreq, leftMinRad)
rightRadSlope, rightRadIntercept, rV, pV, rightErr = linregress(minFreq, rightMinRad)
leftMeanFunc = NP.poly1d(NP.polyfit(minFreq,leftMinRad,2))
fittedLeftMean = leftMeanFunc(minFreq)
rightMeanFunc = NP.poly1d(NP.polyfit(minFreq,rightMinRad,2))
fittedRightMean = rightMeanFunc(minFreq)

fig = PL.figure(figsize = (6,3))
spl = fig.add_subplot(1,1,1)
spl.plot(minFreq,leftMinRad,'.')
spl.plot(minFreq,rightMinRad,'.')
#spl.plot(minFreq,leftRadSlope*NP.array(minFreq) + leftRadIntercept);
#spl.plot(minFreq,rightRadSlope*NP.array(minFreq) + rightRadIntercept);
spl.plot(minFreq,fittedLeftMean)
spl.plot(minFreq,fittedRightMean)
spl.grid()

fig.suptitle('%s %.1f, %.1f arcsec/GHz'%(ephDate, leftRadSlope*60e3, rightRadSlope*60e3))
spl.set_xlabel('MHz')
spl.set_ylabel('arcmin')

del redVis
del redRcp
del redLcp

del redRcpArr
del redLcpArr

#PL.figure(figsize=(8,4))
#PL.imshow(redLcpMean + redRcpMean,aspect=35,origin='lower')
#
#PL.figure(figsize=(8,4))
#for f in NP.linspace(0,31,32,dtype='int'):
#    PL.plot(redLcpMean[f] + redRcpMean[f])

PL.figure()
f193 = fits.open('/home/svlesovoi/sunpy/data/aia20200509_000600_0193.fits')
#f193 = fits.open('/home/svlesovoi/sunpy/data/aia20200510_000300_0193.fits')
x_size = f193[0].header['NAXIS1'] * f193[0].header['CDELT1']
y_size = f193[0].header['NAXIS2'] * f193[0].header['CDELT2']
for f in range(len(minFreq)):
    xx =  x_size/2*NP.sin(leftGp[f] + pAngle)
    yy = -y_size/2*NP.cos(leftGp[f] + pAngle)
    PL.plot([-xx,xx],[-yy,yy],color='green',linewidth=0.3)
    xx =  x_size/2*NP.sin(rightGp[f] + pAngle)
    yy = -y_size/2*NP.cos(rightGp[f] + pAngle)
    PL.plot([-xx,xx],[-yy,yy],color='red',linewidth=0.3)

PL.imshow(f193[0].data,vmax=50,extent=[-x_size/2,x_size/2,y_size/2,-y_size/2])

#del red1RcpMean
#del red1LcpMean
#del red2RcpMean
#del red2LcpMean
