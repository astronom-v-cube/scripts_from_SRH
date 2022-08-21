#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  2 03:37:02 2021

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
from skimage import transform
from skimage import measure
from skimage.feature import canny
from skimage.transform import warp, AffineTransform
from matplotlib.patches import Ellipse, Circle
from zeep import Client

def arcmin_format(xy, pos):
  return '%2d' % ((xy - srh_x_size/2) * resultArcsecPerPixel / 60);

#fileNames = [
#    'SRH36_temp_20210503_1/srh_20210503T013622_clean_image.fit',
#    'SRH36_temp_20210503_1/srh_20210503T013715_clean_image.fit',
#    'SRH36_temp_20210503_1/srh_20210503T013808_clean_image.fit',
#    'SRH36_temp_20210503_1/srh_20210503T013834_clean_image.fit',
#    'SRH36_temp_20210503_1/srh_20210503T013900_clean_image.fit',
#    'SRH36_temp_20210503_1/srh_20210503T013927_clean_image.fit',
#    'SRH36_temp_20210503_1/srh_20210503T013953_clean_image.fit',
#    'SRH36_temp_20210503_1/srh_20210503T014020_clean_image.fit'
#        ]
fileNames = [
    'SRH36_temp_20210509_1/srh_20210509T014748_clean_image.fit',
    'SRH36_temp_20210509_1/srh_20210509T014814_clean_image.fit'
        ]

lcpCenterX = []
lcpCenterY = []
rcpCenterX = []
rcpCenterY = []
lcpECenterX = []
lcpECenterY = []
rcpECenterX = []
rcpECenterY = []

for fileName in fileNames:
    fd = fits.open(fileName)
    
    lcpImage = fd[0].data[0,0]
    rcpImage = fd[0].data[1,0]
    
    srh_x_size = fd[0].header['NAXIS1']
    srh_y_size = fd[0].header['NAXIS2']
    srh_x_delt = fd[0].header['CDELT1'] / 60 * 2
    srh_y_delt = fd[0].header['CDELT2'] / 60 * 2
    
    x0 = int(srh_x_size//2)
    y0 = int(srh_y_size//2)
    dx = int(srh_x_size//5)
    dy = int(srh_y_size//5)
    meanScale = 0.2
    arcsecPerPixel = fd[0].header['CDELT1']*3600
    resultArcsecPerPixel = 4.911/2
    resultScale = arcsecPerPixel / resultArcsecPerPixel
    
    sunCir = skimage.measure.CircleModel()
    sunEll = skimage.measure.EllipseModel()
    
    lcpMean = NP.mean(lcpImage[x0 - dx:x0 + dx, y0 - dy:y0 + dy])
    rcpMean = NP.mean(rcpImage[x0 - dx:x0 + dx, y0 - dy:y0 + dy])
    #lcpImage[x0 - dx:x0 + dx, y0 - dy:y0 + dy] = lcpMean
    #rcpImage[x0 - dx:x0 + dx, y0 - dy:y0 + dy] = rcpMean
    if False:
        contours = measure.find_contours(lcpImage, meanScale*lcpMean)
        contourLength = []
        for n, contour in enumerate(contours):
            contourLength.append(len(contour))
        contourLength = NP.array(contourLength)
        
        maxInd = NP.argmax(contourLength)
        sunCir.estimate(contours[maxInd])
        sunEll.estimate(contours[maxInd])
        lcp_cxc, lcp_cyc, lcp_ca = sunCir.params
        lcp_exc, lcp_eyc, lcp_ea, lcp_eb, lcp_theta = sunEll.params
        cir_patch = Circle((lcp_cyc,lcp_cxc), lcp_ca, edgecolor='yellow', facecolor='none', linewidth=1.)
        
        contours = measure.find_contours(rcpImage, meanScale*rcpMean)
        contourLength = []
        for n, contour in enumerate(contours):
            contourLength.append(len(contour))
        contourLength = NP.array(contourLength)
        
        maxInd = NP.argmax(contourLength)
        sunCir.estimate(contours[maxInd])
        sunEll.estimate(contours[maxInd])
        rcp_cxc, rcp_cyc, rcp_ca = sunCir.params
        rcp_exc, rcp_eyc, rcp_ea, rcp_eb, rcp_theta = sunEll.params
        cir_patch = Circle((rcp_cyc,rcp_cxc), rcp_ca, edgecolor='green', facecolor='none', linewidth=1.)
        
        lcpCenterX.append(lcp_cxc)
        lcpCenterY.append(lcp_cyc)
        rcpCenterX.append(rcp_cxc)
        rcpCenterY.append(rcp_cyc)
        lcpECenterX.append(lcp_exc)
        lcpECenterY.append(lcp_eyc)
        rcpECenterX.append(rcp_exc)
        rcpECenterY.append(rcp_eyc)
        #cir_patch = Circle((y0,x0), lcp_ca, edgecolor='yellow', facecolor='none', linewidth=1.)
        #ell_patch = Ellipse((y0,x0), 2*lcp_ea, 2*lcp_eb, lcp_theta*180./NP.pi, edgecolor='pink', facecolor='none', linewidth=1.)
        shift = AffineTransform(translation=(x0 - lcp_cxc,y0 - lcp_cyc))
        lcpImage = warp(lcpImage,shift.inverse)
        #fig, pl = PL.subplots()
        #pl.set_title('LCP')
        #pl.imshow(lcpImage,origin='lower')
        #pl.add_patch(cir_patch)
        #pl.add_patch(ell_patch)
        
        #cir_patch = Circle((y0,x0), rcp_ca, edgecolor='green', facecolor='none', linewidth=1.)
        #ell_patch = Ellipse((y0,x0), 2*rcp_ea, 2*rcp_eb, rcp_theta*180./NP.pi, edgecolor='pink', facecolor='none', linewidth=1.)
        shift = AffineTransform(translation=(x0 - rcp_cxc,y0 - rcp_cyc))
        rcpImage = warp(rcpImage,shift.inverse)
        #fig, pl = PL.subplots()
        #pl.set_title('RCP')
        #pl.imshow(rcpImage,origin='lower')
        #pl.add_patch(cir_patch)
        #pl.add_patch(ell_patch)
        
        pAngle = 0
        try:
            client = Client('http://ephemeris.rao.istp.ac.ru/?wsdl')
            result = client.service.Ephemeride('SSRT','sun',fd[0].header['DATE-OBS'])
            pAngle = NP.deg2rad(float(result[0]['PAngle']))
        except:
            pass
        
        scale = AffineTransform(scale=(-resultScale,-resultScale))
        shift = AffineTransform(translation=(-srh_y_size/2,-srh_y_size/2))
        rotate = AffineTransform(rotation = -pAngle)
        back_shift = AffineTransform(translation=(srh_y_size/2,srh_y_size/2))
        
        lcpImage = warp(lcpImage,(shift + (rotate + back_shift)).inverse)
        rcpImage = warp(rcpImage,(shift + (rotate + back_shift)).inverse)
        lcpImage = warp(lcpImage,(shift + (scale + back_shift)).inverse)
        rcpImage = warp(rcpImage,(shift + (scale + back_shift)).inverse)
        
        fig, pl = PL.subplots(figsize=(8,8))
        fig.tight_layout()
        cir_patch = Circle((y0,x0), rcp_ca*resultScale, edgecolor='red', facecolor='none', linewidth=.3)
        
        pl.xaxis.set_major_locator(PL.MultipleLocator(256));
        pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
        pl.xaxis.set_minor_locator(PL.MultipleLocator(64));
        pl.yaxis.set_major_locator(PL.MultipleLocator(256));
        pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
        pl.yaxis.set_minor_locator(PL.MultipleLocator(64));
        pl.set_title('SRH V, %s' % fd[0].header['DATE-OBS'])
        pl.set_xlabel('arcmin')
        pl.set_ylabel('arcmin')
        pl.add_patch(cir_patch)
        pl.imshow(lcpImage - rcpImage,origin='lower',cmap='Greys')
        fig.savefig('srh_V_%s.png'%fd[0].header['DATE-OBS'])
        
        fig, pl = PL.subplots(figsize=(8,8))
        fig.tight_layout()
        cir_patch = Circle((y0,x0), rcp_ca*resultScale, edgecolor='white', facecolor='none', linewidth=.3)
        pl.xaxis.set_major_locator(PL.MultipleLocator(256));
        pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
        pl.xaxis.set_minor_locator(PL.MultipleLocator(64));
        pl.yaxis.set_major_locator(PL.MultipleLocator(256));
        pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
        pl.yaxis.set_minor_locator(PL.MultipleLocator(64));
        pl.set_title('SRH I %s' % fd[0].header['DATE-OBS'])
        pl.set_xlabel('arcmin')
        pl.set_ylabel('arcmin')
        #pl.add_patch(cir_patch)
        pl.imshow(lcpImage + rcpImage,origin='lower',cmap='hot')
        fig.savefig('srh_I_%s.png'%fd[0].header['DATE-OBS'])
    
    stdV = NP.zeros((10,10))
    dL = 1.2
    for i in range(10):
        for j in range(10):
            shift = AffineTransform(translation=((i-5)*dL,(j-5)*dL))
            stdV[i,j]= ((warp(rcpImage,shift.inverse) - lcpImage)**2).mean()
    
    rlShift = NP.where(stdV==stdV.min())
    rlShift = (rlShift[0]%lcpImage.shape[0], rlShift[1]%lcpImage.shape[1])
    rlShift = [(rlShift[0][0]-5)*dL,(rlShift[1][0]-5)*dL]
    PL.figure()
    shift = AffineTransform(translation=rlShift)
    rcpImage = warp(rcpImage,shift.inverse)
#    PL.imshow(lcpImage - rcpImage,vmin=-0.0001,vmax=0.0001,origin='lower',cmap='Greys')
    PL.imshow(rcpImage - lcpImage,vmin=-0.0001,vmax=0.0001,origin='lower',cmap='Greys')

    fd.close()
    