#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 09:59:07 2020

@author: mariagloba
"""

import numpy as NP
from srhFitsFile_doubleBase import SrhFitsFile
from skimage.transform import warp, AffineTransform
from astropy.io import fits
from astropy.time import Time, TimeDelta
import base2uvw as bl2uvw
import pylab as PL
from scipy.signal import argrelextrema
from scipy.optimize import minimize
from sunpy import coordinates
import os
from skimage.transform import hough_circle, hough_circle_peaks, hough_ellipse
from skimage.feature import canny
import skimage
from matplotlib.patches import Ellipse
import time

class srhImage():
    def buildSPhase(self):
        self.sLcpPhaseCorrection[:] = 0.
        self.sRcpPhaseCorrection[:] = 0.
        for j in range(16):
            for i in range(16):
                self.sLcpPhaseCorrection[i] +=  NP.deg2rad(self.sPhaseCoefsLcp[j] * (-1)**(i // (j + 1))) 
                self.sRcpPhaseCorrection[i] +=  NP.deg2rad(self.sPhaseCoefsRcp[j] * (-1)**(i // (j + 1))) 
        self.srhFits.changeSouthPhase(self.sLcpPhaseCorrection, self.sRcpPhaseCorrection)
     
    def lm2heliocentric(self):
        scaling = self.srhFits.getPQScale(self.uvSize, NP.deg2rad(self.arcsecPerPixel * (self.uvSize - 1)/3600.))
        scale = AffineTransform(scale=(self.uvSize/scaling[0], self.uvSize/scaling[1]))
        shift = AffineTransform(translation=(-self.uvSize/2,-self.uvSize/2))
        rotate = AffineTransform(rotation = self.pAngle)
        matrix = AffineTransform(matrix = self.srhFits.getPQ2HDMatrix())
        back_shift = AffineTransform(translation=(self.uvSize/2,self.uvSize/2))
    
        dataResult0 = warp(self.srhFits.lcp.real,(shift + (scale + back_shift)).inverse)
        lcpData = warp(dataResult0,(shift + (matrix + back_shift)).inverse)
        lcpData = warp(lcpData,(shift + (rotate + back_shift)).inverse)
        dataResult0 = warp(self.srhFits.rcp.real,(shift + (scale + back_shift)).inverse)
        rcpData = warp(dataResult0,(shift + (matrix + back_shift)).inverse)
        rcpData = warp(rcpData,(shift + (rotate + back_shift)).inverse)
        
        self.imagesLcp = lcpData
        self.imagesRcp = rcpData
    
    def buildRawImage(self):
        print('making raw images...')
        self.srhFits.averageCalib = True
        self.srhFits.updateAntennaPhase()
        self.srhFits.updateAntennaAmplitude()
        self.srhFits.vis2uv(0, phaseCorrect=True, amplitudeCorrect=True, PSF=False, average = self.srhFits.dataLength)
        self.srhFits.uv2lmImage()
        self.imagesLcpLm = self.srhFits.lcp.real.copy()
        self.imagesRcpLm = self.srhFits.rcp.real.copy()
        
        self.lm2heliocentric()
        
    def centeringCoarse(self):
        print('centering...')
        PL.ioff()
        edges = canny(self.imagesLcp, sigma=3, low_threshold=0, high_threshold=2e-4)
        hough_radii = NP.arange(90, 120, 2)
        hough_res = hough_circle(edges, hough_radii)
        accums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii,
                                               total_num_peaks=1)
        self.cxcyLcp = NP.array([cx[0], cy[0]])
        
        edges = canny(self.imagesRcp, sigma=3, low_threshold=0, high_threshold=2e-4)
        hough_radii = NP.arange(90, 120, 2)
        hough_res = hough_circle(edges, hough_radii)
        accums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii,
                                               total_num_peaks=1)
        self.cxcyRcp = NP.array([cx[0], cy[0]])
#            fig, ax = PL.subplots(ncols=1, nrows=1, figsize=(10, 10))
#            
#            for center_y, center_x, radius in zip(cy, cx, radii):
#    
#                ax.imshow(self.imagesLcp[self.freq])
#                c = PL.Circle((center_x, center_y), radius, color='red', linewidth=2, fill=False)
#                ax.add_patch(c)
#            PL.title(str(self.srhFits.freqList[self.freq]))
#            PL.savefig('/home/mariagloba/Work/fits/20200815/centering_test/' + str(self.srhFits.freqList[freq]) + '.png')
#            PL.close()
            
    def centeringFine_circle(self):
        print('centering...')
        PL.ioff()
        size = self.uvSize//2
        d = 70
        O = size//2
        cxl, cyl = self.cxcyLcp.astype(int)
        data = self.imagesLcp[cyl - size//2:cyl + size//2, cxl - size//2:cxl + size//2]
        mean = NP.mean(data[O - d:O + d, O - d:O + d])
        edges = canny(NP.clip(data,mean/4,mean/4 + mean/40), sigma=3)
        hough_radii = NP.arange(90, 120, 1)
        hough_res = hough_circle(edges, hough_radii)
        accums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii,
                                               total_num_peaks=1)
        
        self.cxcyLcp[0] += cy[0] - O
        self.cxcyLcp[1] += cx[0] - O
        
        transVector = AffineTransform(translation = -self.cxcyLcp + 256)
        self.imagesLcp = warp(self.imagesLcp, transVector.inverse)
        self.clippedImagesLcp = self.imagesLcp[size-O:size+O, size-O:size+O]

       
        cxr, cyr = self.cxcyRcp.astype(int)
        data = self.imagesRcp[cyr - size//2:cyr + size//2, cxr - size//2:cxr + size//2]
        mean = NP.mean(data[O - d:O + d, O - d:O + d])
        edges = canny(NP.clip(data,mean/4,mean/4 + mean/40), sigma=3)
        hough_radii = NP.arange(90, 120, 1)
        hough_res = hough_circle(edges, hough_radii)
        accums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii,
                                               total_num_peaks=1)
        
        self.cxcyRcp[0] += cy[0] - O
        self.cxcyRcp[1] += cx[0] - O
        
        transVector = AffineTransform(translation = -self.cxcyRcp + 256)
        self.imagesRcp = warp(self.imagesRcp, transVector.inverse)
        self.clippedImagesRcp = self.imagesRcp[size-O:size+O, size-O:size+O]
        
#        PL.clf()
#        PL.imshow(self.clippedImagesLcp+self.clippedImagesRcp, cmap = 'ocean', vmin = -84000, vmax = 84000)
#        PL.title(str(self.srhFits.freqList[self.freq]))
#        PL.savefig('/home/mariagloba/Work/fits/20200805/centering_test_circle/final_I_' + str(self.srhFits.freqList[self.freq]) + '.png')
#       
#        PL.clf()
#        PL.imshow(self.clippedImagesLcp-self.clippedImagesRcp, cmap = 'ocean')
#        PL.title(str(self.srhFits.freqList[self.freq]))
#        PL.savefig('/home/mariagloba/Work/fits/20200805/centering_test_circle/final_V_' + str(self.srhFits.freqList[self.freq]) + '.png')
#           
        
#            fig, ax = PL.subplots(ncols=1, nrows=1, figsize=(10, 10))
#            
#            for center_y, center_x, radius in zip(cy, cx, radii):
#    
#                ax.imshow(self.imagesLcp[freq])
#                c = PL.Circle((center_x, center_y), radius, color='red', linewidth=2, fill=False)
#                ax.add_patch(c)
#            PL.title(str(self.srhFits.freqList[freq]))
#            PL.savefig('/home/mariagloba/Work/fits/20200815/centering_test/' + str(self.srhFits.freqList[freq]) + '.png')
#            PL.close()
        

            
    def findPhase(self):
        print('calculating phase stair...')
#        PL.ioff()
#        PL.figure()
        d = 70
        self.lcpMaxTrace = []
        self.rcpMaxTrace = []
#            self.srhFits.setSizeOfUv(256)
        cxl, cyl = self.cxcyLcp.astype(int)
        cxr, cyr = self.cxcyRcp.astype(int)
        for sPhaseInd in range(-18,18):
            self.sLcpPhaseCorrection[:] = NP.deg2rad(sPhaseInd*10)
            self.sRcpPhaseCorrection[:] = NP.deg2rad(sPhaseInd*10)
            self.srhFits.changeSouthPhase(self.sLcpPhaseCorrection, self.sRcpPhaseCorrection)
            self.srhFits.vis2uv(0, phaseCorrect=True, amplitudeCorrect=True, PSF=False, average = self.srhFits.dataLength);
            self.srhFits.uv2lmImage()
            self.lm2heliocentric()
            self.lcpMaxTrace.append(self.imagesLcp[cyl-d:cyl+d, cxl-d:cxl+d].mean())
            self.rcpMaxTrace.append(self.imagesRcp[cyr-d:cyr+d, cxr-d:cxr+d].mean())
#            self.srhFits.setSizeOfUv(self.uvSize)
        
        phaseIndLcp = int(10*(NP.argmax(self.lcpMaxTrace) - 18) + .5)
        phaseIndRcp = int(10*(NP.argmax(self.rcpMaxTrace) - 18) + .5)
      
        self.sPhaseCoefsLcp[15] = phaseIndLcp
        self.sPhaseCoefsRcp[15] = phaseIndRcp
        
        self.lcpMaxTrace = []
        self.rcpMaxTrace = []
        
        for sPhaseInd in range(-10,10):
            self.lInd = (phaseIndLcp + sPhaseInd + 180) % 360 - 180
            self.rInd = (phaseIndRcp + sPhaseInd + 180) % 360 - 180
            self.sLcpPhaseCorrection[:] = NP.deg2rad(self.lInd)
            self.sRcpPhaseCorrection[:] = NP.deg2rad(self.rInd)
            self.srhFits.changeSouthPhase(self.sLcpPhaseCorrection, self.sRcpPhaseCorrection)
            self.srhFits.vis2uv(0, phaseCorrect=True, amplitudeCorrect=True, PSF=False, average = self.srhFits.dataLength);
            self.srhFits.uv2lmImage()
            self.hLcp, self.bLcp = NP.histogram(self.srhFits.lcp.real, bins = 101)
            self.hRcp, self.bRcp = NP.histogram(self.srhFits.rcp.real, bins = 101)
            self.lcpMaxTrace.append(NP.max(self.hLcp))
            self.rcpMaxTrace.append(NP.max(self.hLcp))
        
        phaseIndLcpFine = int(NP.argmax(self.lcpMaxTrace) - 10 + phaseIndLcp)
        phaseIndRcpFine = int(NP.argmax(self.rcpMaxTrace) - 10 + phaseIndRcp)
        self.sPhaseCoefsLcp[15] = phaseIndLcpFine
        self.sPhaseCoefsRcp[15] = phaseIndRcpFine
            
        self.buildSPhase()
        self.srhFits.vis2uv(0, phaseCorrect=True, amplitudeCorrect=True, PSF=False, average = self.srhFits.dataLength)
        self.srhFits.uv2lmImage()
        self.imagesLcpLm = self.srhFits.lcp.real.copy()
        self.imagesRcpLm = self.srhFits.rcp.real.copy()
        self.lm2heliocentric()
        
#            PL.clf()
#            PL.imshow(self.imagesLcp[freq][cyl-125:cyl+125, cxl-125:cxl+125])
#            PL.title(str(self.srhFits.freqList[freq]))
#            PL.savefig('/home/mariagloba/Work/fits/20200815/findPhase_test/' + str(self.srhFits.freqList[freq]) + '_lcp.png')
#            PL.clf()
#            PL.imshow(self.imagesRcp[freq][cyr-125:cyr+125, cxr-125:cxr+125])
#            PL.title(str(self.srhFits.freqList[freq]))
#            PL.savefig('/home/mariagloba/Work/fits/20200815/findPhase_test/' + str(self.srhFits.freqList[freq]) + '_rcp.png')
#        PL.close()
        
    def calibrateBrightness(self):
        print('calibrating brightness...')
        h_l, b_l = NP.histogram(self.imagesLcpLm, bins = 101)
        h_r, b_r = NP.histogram(self.imagesRcpLm, bins = 101)
        sun_l = b_l[NP.argmax(h_l[60:])+60]
        sun_r = b_r[NP.argmax(h_r[60:])+60]
        self.calCoefLcp = 16000./sun_l
        self.calCoefRcp = 16000./sun_r
        self.imagesLcp *= self.calCoefLcp
        self.imagesRcp *= self.calCoefRcp

    def centeringFine(self):
        print('centering...')
        size = self.uvSize//2
        d = 70
        O = size//2
        sunEll = skimage.measure.EllipseModel()
        cxl, cyl = self.cxcyLcp.astype(int)
        data = self.imagesLcp[cyl - size//2:cyl + size//2, cxl - size//2:cxl + size//2]
        mean = NP.mean(data[O - d:O + d, O - d:O + d])
        edge = skimage.filters.sobel(NP.clip(data,mean/4,mean/4 + mean/40))
        ellInd = NP.where(edge > .5*edge.max())
        ellXY = NP.zeros((len(ellInd[1]),2))
        ellXY[:,0] = ellInd[0]
        ellXY[:,1] = ellInd[1]
        ellXY[:,0] = (ellXY[:,0] - O)# * srhFits[0].header['CDELT1'] / 60 * 2
        ellXY[:,1] = (ellXY[:,1] - O)# * srhFits[0].header['CDELT2'] / 60 * 2
        sunEll.estimate(ellXY)
        self.xc_lcp, self.yc_lcp, b, a, theta = sunEll.params
        self.cxcyLcp[0] += self.yc_lcp
        self.cxcyLcp[1] += self.xc_lcp
        
        transVector = AffineTransform(translation = -self.cxcyLcp + 256)
        self.imagesLcp = warp(self.imagesLcp, transVector.inverse)
        self.clippedImagesLcp = self.imagesLcp[size-O:size+O, size-O:size+O]

            
#            PL.figure()
#            fig, pl = PL.subplots()
#            pl.imshow(data)
#            ell_patch = Ellipse((O+self.yc_lcp,O+self.xc_lcp), 2*a, 2*b, -theta*180/NP.pi, edgecolor='red', facecolor='none', linewidth=1.)
#            pl.add_patch(ell_patch)
#            pl.grid()
#            PL.title(str(self.srhFits.freqList[freq]))
#            PL.savefig('/home/mariagloba/Work/fits/20200815/centering_test/lcp_' + str(self.srhFits.freqList[freq]) + '.png')
#            PL.close()
            
        
        cxr, cyr = self.cxcyRcp.astype(int)
        data = self.imagesRcp[cyr - size//2:cyr + size//2, cxr - size//2:cxr + size//2]
        mean = NP.mean(data[O - d:O + d, O - d:O + d])
        edge = skimage.filters.sobel(NP.clip(data,mean/4,mean/4 + mean/40))
        ellInd = NP.where(edge > .5*edge.max())
        ellXY = NP.zeros((len(ellInd[1]),2))
        ellXY[:,0] = ellInd[0]
        ellXY[:,1] = ellInd[1]
        ellXY[:,0] = (ellXY[:,0] - O)# * srhFits[0].header['CDELT1'] / 60 * 2
        ellXY[:,1] = (ellXY[:,1] - O)# * srhFits[0].header['CDELT2'] / 60 * 2
        sunEll.estimate(ellXY)
        self.xc_rcp, self.yc_rcp, b, a, theta = sunEll.params
        self.cxcyRcp[0] += self.yc_rcp
        self.cxcyRcp[1] += self.xc_rcp
        
        transVector = AffineTransform(translation = -self.cxcyRcp + 256)
        self.imagesRcp = warp(self.imagesRcp, transVector.inverse)
        
        self.vCorrection()
        
        self.clippedImagesRcp = self.imagesRcp[size-O:size+O, size-O:size+O]
#
#            PL.clf()
#            PL.imshow(self.clippedImagesLcp[self.freq]+self.clippedImagesRcp[self.freq], cmap = 'ocean', vmin = -84000, vmax = 84000)
#            PL.title(str(self.srhFits.freqList[self.freq]))
#            PL.savefig('/home/mariagloba/Work/fits/20200805/centering_test_mini/final_I_' + str(self.srhFits.freqList[freq]) + '.png')
#
#            PL.clf()
#            PL.imshow(self.clippedImagesLcp[freq] - self.clippedImagesRcp[freq], cmap = 'ocean')#, vmin = -42000, vmax = 42000)
#            PL.title(str(self.srhFits.freqList[freq]))
#            PL.savefig('/home/mariagloba/Work/fits/20200805/centering_test_mini/final_V_' + str(self.srhFits.freqList[freq]) + '.png')

#
#            fig, pl = PL.subplots()
#            pl.imshow(data)
#            ell_patch = Ellipse((O+self.yc_rcp,O+self.xc_rcp), 2*a, 2*b, -theta*180/NP.pi, edgecolor='red', facecolor='none', linewidth=1.)
#            pl.add_patch(ell_patch)
#            pl.grid()
#            PL.title(str(self.srhFits.freqList[freq]))
#            PL.savefig('/home/mariagloba/Work/fits/20200815/centering_test/rcp_' + str(self.srhFits.freqList[freq]) + '.png')
#            PL.close()
#        PL.close()
        
    def mean_V(self, shift):
        size = self.uvSize//2
        O = size//2
        transVector = AffineTransform(translation = shift)
        rcp = warp(self.imagesRcp, transVector.inverse)
        return NP.mean(NP.abs((self.imagesLcp-rcp)[size-O:size+O, size-O:size+O]))
        
    def vCorrection(self):
        print('correcting lcp-rcp shift...')
        self.res = minimize(self.mean_V, [0,0])
        transVector = AffineTransform(translation = self.res.x)
        self.imagesRcp = warp(self.imagesRcp, transVector.inverse)

#    def clipImages(self):
#        
    def saveAsFits(self):
        fileDate = self.srhFits.hduList[0].header['DATE-OBS'].replace('/', '')
        fileTime = self.srhFits.hduList[0].header['TIME-OBS'].replace(':','').split('.')[0]
        fitsNameLcp = fileDate + '_' + fileTime + '_' + str(int(self.srhFits.freqList[self.freq])) + '_LCP.fit'
        fitsNameRcp = fileDate + '_' + fileTime + '_' + str(int(self.srhFits.freqList[self.freq])) + '_RCP.fit'
        print(fitsNameLcp)
        pHeader = fits.Header();
        pHeader['DATE-OBS']     = self.srhFits.hduList[0].header['DATE-OBS'].replace('/', '-')
        pHeader['T-OBS']        = self.srhFits.hduList[0].header['TIME-OBS']
        pHeader['INSTRUME']     = self.srhFits.hduList[0].header['INSTRUME']
        pHeader['ORIGIN']       = self.srhFits.hduList[0].header['ORIGIN']
        pHeader['OBS-LAT']      = self.srhFits.hduList[0].header['OBS-LAT']
        pHeader['OBS-LONG']     = self.srhFits.hduList[0].header['OBS-LONG']
        pHeader['OBS-ALT']      = self.srhFits.hduList[0].header['OBS-ALT']
        pHeader['FR_CHAN']      = self.srhFits.hduList[0].header['FR_CHAN']
        pHeader['FREQUENC']     = ('%3.3f') % (self.srhFits.freqList[self.freq]*1e6)
        pHeader['CDELT1']       = self.arcsecPerPixel
        pHeader['CDELT2']       = self.arcsecPerPixel
        pHeader['CRPIX1']       = self.uvSize // 2
        pHeader['CRPIX2']       = self.uvSize // 2
        pHeader['CTYPE1']       = 'HPLN-TAN'
        pHeader['CTYPE2']       = 'HPLT-TAN'
        pHeader['CUNIT1']       = 'arcsec'
        pHeader['CUNIT2']       = 'arcsec'
        
        
        pHdu = fits.PrimaryHDU(NP.flip(self.clippedImagesLcp, 0), header=pHeader);
        hduList = fits.HDUList([pHdu]);
        hduList.writeto(os.path.join(self.outPath, fitsNameLcp), clobber=True);
        hduList.close();
        
        pHdu = fits.PrimaryHDU(NP.flip(self.clippedImagesRcp, 0), header=pHeader);
        hduList = fits.HDUList([pHdu]);
        hduList.writeto(os.path.join(self.outPath, fitsNameRcp), clobber=True);
        
        hduList.close();
    
    def __init__(self, filename, freq):
        self.uvSize = 512
        self.mfFitsName = filename
        self.freq = freq
        self.imagesLcpLm = NP.zeros((512, 512))
        self.imagesRcpLm = NP.zeros((512, 512))
        self.imagesLcp = NP.zeros((512, 512))
        self.imagesRcp = NP.zeros((512, 512))
        self.clippedImagesLcp = NP.zeros((256, 256))
        self.clippedImagesRcp = NP.zeros((256, 256))
        self.ewLcpPhaseCorrection = NP.zeros(32)
        self.ewRcpPhaseCorrection = NP.zeros(32)
        self.sLcpPhaseCorrection = NP.zeros(16)
        self.sRcpPhaseCorrection = NP.zeros(16)
        self.pAngle = 0
        self.arcsecPerPixel = 4.91104 * 2
        self.cxcyLcp = NP.zeros(2)
        self.cxcyRcp = NP.zeros(2)
        self.outPath = './'
        
        self.srhFits = SrhFitsFile(filename, self.uvSize)
        self.pAngle = NP.deg2rad(coordinates.get_sun_P(self.srhFits.dateObs).to_value())
        self.srhFits.getHourAngle(self.srhFits.dataLength//2)
        self.srhFits.setFrequencyChannel(self.freq)
        
        self.ewPhaseCoefsLcp = NP.zeros(16)
        self.ewPhaseCoefsRcp = NP.zeros(16)
        self.sPhaseCoefsLcp = NP.zeros(16)
        self.sPhaseCoefsRcp = NP.zeros(16)
        
    def makeImages(self):
        self.buildRawImage()
        self.centeringCoarse()
        self.findPhase()
        self.calibrateBrightness()
        self.centeringFine()
#        self.vCorrection()
        self.saveAsFits()
        PL.ion()
        