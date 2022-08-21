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
import scipy.signal
import scipy.constants

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

#fileNames = [
#    'SRH36_temp_20210509_1/srh_20210509T014748_clean_image.fit',
#    'SRH36_temp_20210509_1/srh_20210509T014814_clean_image.fit',
#    'SRH36_temp_20210509_1/srh_20210509T014841_clean_image.fit',
#    'SRH36_temp_20210509_1/srh_20210509T014907_clean_image.fit',
#    'SRH36_temp_20210509_1/srh_20210509T014933_clean_image.fit'
#        ]
#
#fileNames = [
#    'SRH36_temp_20210512_1/srh_20210512T055454_clean_image.fit'
#        ]
#fileNames = [
#    'SRH36_temp_20210512_2/srh_20210512T014516_clean_image.fit'
#        ]
#fileNames = [
#    'SRH36_temp_20210513_1/srh_20210513T034149_clean_image.fit'
#        ]
#fileNames = [
#    'SRH36_temp_20210513_2/srh_20210513T072116_clean_image.fit',
#    'SRH36_temp_20210513_2/srh_20210513T072142_clean_image.fit',
#    'SRH36_temp_20210513_2/srh_20210513T072209_clean_image.fit',
#    'SRH36_temp_20210513_2/srh_20210513T072235_clean_image.fit'
#        ]
#
#fileNames = [
#    'SRH36_temp_20210511_1/srh_20210511T024852_clean_image.fit'
#        ]

#fileNames = [
#    'SRH36_temp_20210510_1/srh_20210510T024749_clean_image.fit'
#        ]
#fileNames = [
#    'SRH36_temp_20210516_1/srh_20210516T031151_clean_image.fit'
#        ]
#fileNames = [
#    'SRH36_temp_20210515_1/srh_20210515T031232_clean_image.fit'
#        ]
#fileNames = [
#    'SRH36_temp_20210520_1/srh_20210520T030715_clean_image.fit'
#        ]
#psfName = 'SRH36_temp_20210520_1/srh_20210520T030715_clean_psf.fit'
#
fileNames = [
    'SRH36_temp_20210522_1/srh_20210522T015224_clean_image.fit'
        ]
psfName = 'SRH36_temp_20210522_1/srh_20210522T015224_clean_psf.fit'

filter_dL = 21
filter_arg = NP.linspace(-1,1,filter_dL)
filterX, filterY = NP.meshgrid(filter_arg, filter_arg)
imageFilter = NP.exp(-((filterX/5.7)**2 + (filterY/5.7)**2))

iImageMaxX = []
iImageMaxY = []
flux107 = []
#Tb_2800 = 27100
Tb_2800 = 27100

psf_fd = fits.open(psfName)
psfData = psf_fd[0].data[0][0]

for fileName in fileNames:
    fd = fits.open(fileName)
    
    pAngle = 0
    try:
        client = Client('http://ephemeris.rao.istp.ac.ru/?wsdl')
        result = client.service.Ephemeride('SSRT','sun',fd[0].header['DATE-OBS'])
        pAngle = NP.deg2rad(float(result[0]['PAngle']))
    except:
        pass
    
    lcpImage = fd[0].data[0,0]
    rcpImage = fd[0].data[1,0]
    
    sunMask = NP.zeros_like(lcpImage)
    outRadius = 400
    inRadius = 290
    sunX0 = lcpImage.shape[0]/2
    sunY0 = lcpImage.shape[0]/2
    for i in range(lcpImage.shape[0]):
        for j in range(lcpImage.shape[0]):
            curRadius = NP.sqrt((i - sunX0)**2 + (j - sunX0)**2)
            if curRadius < outRadius:
                sunMask[i,j] = 1.

    srh_x_size = fd[0].header['NAXIS1']
    srh_y_size = fd[0].header['NAXIS2']
    srh_x_delt = fd[0].header['CDELT1'] / 60 * 2
    srh_y_delt = fd[0].header['CDELT2'] / 60 * 2
    
    x0 = int(srh_x_size//2)
    y0 = int(srh_y_size//2)
    dx = int(srh_x_size//5)
    dy = int(srh_y_size//5)
    meanScale = 0.15
    arcsecPerPixel = fd[0].header['CDELT1']*3600
    resultArcsecPerPixel = 4.911/2
    resultScale = arcsecPerPixel / resultArcsecPerPixel
    
    sunCir = skimage.measure.CircleModel()
    sunEll = skimage.measure.EllipseModel()
    
    lcpMean = NP.mean(lcpImage[x0 - dx:x0 + dx, y0 - dy:y0 + dy])
    rcpMean = NP.mean(rcpImage[x0 - dx:x0 + dx, y0 - dy:y0 + dy])
    stdV = NP.zeros((10,10))
    dL = 1.2
    for i in range(10):
        for j in range(10):
            shift = AffineTransform(translation=((i-5)*dL,(j-5)*dL))
            stdV[i,j]= ((warp(rcpImage,shift.inverse) - lcpImage)**2).mean()
    
    rlShift = NP.where(stdV==stdV.min())
    rlShift = (rlShift[0]%lcpImage.shape[0], rlShift[1]%lcpImage.shape[1])
    rlShift = [(rlShift[0][0]-5)*dL,(rlShift[1][0]-5)*dL]
    shift = AffineTransform(translation=rlShift)
    rcpImage = warp(rcpImage,shift.inverse)
    iImage = rcpImage + lcpImage
    vImage = rcpImage - lcpImage
    
    meanI = (iImage*sunMask).mean()
#    iImageHist = NP.histogram(NP.clip(iImage*sunMask,-meanI,5*meanI),bins=1000)
    iImageHist = NP.histogram(iImage,bins=1000)
    maxInds = scipy.signal.find_peaks(iImageHist[0],prominence=2000)
    iImage -= iImageHist[1][maxInds[0][0]]
    iImage /= (iImageHist[1][maxInds[0][1]] - iImageHist[1][maxInds[0][0]])
    iImage *= Tb_2800
    vImage /= (iImageHist[1][maxInds[0][1]] - iImageHist[1][maxInds[0][0]])
    vImage *= Tb_2800
    
    iMean = NP.mean(iImage[x0 - dx:x0 + dx, y0 - dy:y0 + dy])
    if True:
        scale = AffineTransform(scale=(-resultScale,-resultScale))
        shift = AffineTransform(translation=(-srh_y_size/2,-srh_y_size/2))
        rotate = AffineTransform(rotation = -pAngle)
        back_shift = AffineTransform(translation=(srh_y_size/2,srh_y_size/2))
        
        iImage = warp(iImage,(shift + (rotate + back_shift)).inverse)
        vImage = warp(vImage,(shift + (rotate + back_shift)).inverse)
        iImage = warp(iImage,(shift + (scale + back_shift)).inverse)
        vImage = warp(vImage,(shift + (scale + back_shift)).inverse)
        
        conImage = scipy.signal.fftconvolve(iImage,imageFilter) / filter_dL**2
        conImage = conImage[filter_dL//2:filter_dL//2 + srh_x_size,filter_dL//2:filter_dL//2 + srh_y_size]
#        contours = measure.find_contours(conImage, meanScale*iMean)
        contours = measure.find_contours(conImage, 17000)
        contourLength = []
        for n, contour in enumerate(contours):
            contourLength.append(len(contour))
        contourLength = NP.array(contourLength)
        maxInd = NP.argmax(contourLength)
        sunCir.estimate(contours[maxInd])
        i_cx0, i_cy0, i_cR = sunCir.params

        iShift = AffineTransform(translation=(y0 + .5 - i_cy0,x0 + .5 - i_cx0))
        iImage = warp(iImage,iShift.inverse)
        vImage = warp(vImage,iShift.inverse)

        fig, pl = PL.subplots(figsize=(8,8))
        fig.tight_layout()
        cir_patch = Circle((y0,x0), i_cR, edgecolor='red', facecolor='none', linewidth=1.3)
        iCir_patch = Circle((i_cy0,i_cx0), i_cR, edgecolor='green', facecolor='none', linewidth=1.3)
        pl.xaxis.set_major_locator(PL.MultipleLocator(256));
        pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
        pl.xaxis.set_minor_locator(PL.MultipleLocator(64));
        pl.yaxis.set_major_locator(PL.MultipleLocator(256));
        pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
        pl.yaxis.set_minor_locator(PL.MultipleLocator(64));
        pl.text(10,25,'SRH V, %s' % fd[0].header['DATE-OBS'],color='white',fontsize='14')
        pl.set_xlabel('arcmin')
        pl.set_ylabel('arcmin')
        pl.imshow(vImage,cmap='Greys',vmin=-0.3e6,vmax=0.3e6,origin='lower')
        pl.add_patch(cir_patch)
        fig.savefig('srh_V_%s.png'%fd[0].header['DATE-OBS'])
        
        maxI = NP.where(iImage==iImage.max())
        iImageMaxX.append(maxI[0][0]%iImage.shape[0])
        iImageMaxY.append(maxI[1][0]%iImage.shape[1])
#        flux107.append(2*(iImage*sunMask).mean()*scipy.constants.k/(0.107**2)*((1024*NP.deg2rad(fd[0].header['CDELT1']))**2)/1e-22)
        flux107.append(2*(iImage*sunMask).sum()*scipy.constants.k/(0.107**2)*(NP.deg2rad(fd[0].header['CDELT1'])**2)/1e-22)

        fig, pl = PL.subplots(figsize=(8,8))
        fig.tight_layout()
        cir_patch = Circle((y0,x0), i_cR, edgecolor='red', facecolor='none', linewidth=.3)
        pl.xaxis.set_major_locator(PL.MultipleLocator(256));
        pl.xaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
        pl.xaxis.set_minor_locator(PL.MultipleLocator(64));
        pl.yaxis.set_major_locator(PL.MultipleLocator(256));
        pl.yaxis.set_major_formatter(PL.FuncFormatter(arcmin_format));
        pl.yaxis.set_minor_locator(PL.MultipleLocator(64));
        pl.text(10,25,'SRH I %s, S=%.1f s.f.u.' % (fd[0].header['DATE-OBS'], flux107[-1]),color='white',fontsize='14')
        pl.set_xlabel('arcmin')
        pl.set_ylabel('arcmin')
#        pl.add_patch(cir_patch)
        pl.imshow(iImage,cmap='hot',vmin=0.,vmax=3e5,origin='lower')
#        pl.contour(vImage,cmap='bwr',levels=[-20e4,-10e4,-5e4,-2e4,-1e4,1e4,2e4,5e4,10e4,20e4],linewidths=0.5)
        pl.contour(NP.roll(NP.roll(psfData,-srh_x_size//2-100,axis=0),srh_y_size//2+100,axis=1),levels=[0.5],cmap='gray_r')
    
        fig.savefig('srh_I_%s.png'%fd[0].header['DATE-OBS'])

    fd.close()
    