#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 14:24:59 2022

@author: maria
"""
print('Start')
from casatasks import tclean, exportfits, flagdata, rmtables
from casatools import image as IA
from casacore.images import image
import numpy as NP
import os, fnmatch
import scipy.signal
from srhFitsFile612 import SrhFitsFile
import srh612MS2
import skimage.measure
from astropy.io import fits
from skimage.transform import warp, AffineTransform
import pylab as PL
import datetime

def ihhmm_format(t):
  hh = int(t / 3600.)
  t -= hh*3600.
  mm = int(t / 60.)
  t -= mm*60.
  ss = int(t)
  return '%02d:%02d:%02d' % (hh,mm,ss)

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;


def createDiskMask(N, cell, arcsecRadius):
    radius = int(arcsecRadius/cell +0.5)
    model = NP.zeros((N,N))
    for i in range(N):
        for j in range(N):
            x=i - N/2
            y=j - N/2
            if (NP.sqrt(x**2 + y**2) < radius):
                model[i, j] = 1
    return model

def saveFitsImages(iImage, vImage, restoring_beam, saveCleanFitsPath, srhRawFits, frequencyIndex, scan):
    resultArcsecPerPixel = 4.9
    srh_x_size = 1648/2
    
    a,b,ang = restoring_beam['major']['value'],restoring_beam['minor']['value'],restoring_beam['positionangle']['value']
    fitsTime = ihhmm_format(srhRawFits.freqTime[frequencyIndex, scan])

    pHeader = fits.Header();
    pHeader['DATE-OBS']     = srhRawFits.hduList[0].header['DATE-OBS']
    pHeader['T-OBS']        = fitsTime#srhRawFits.hduList[0].header['TIME-OBS']
    pHeader['INSTRUME']     = srhRawFits.hduList[0].header['INSTRUME']
    pHeader['ORIGIN']       = srhRawFits.hduList[0].header['ORIGIN']
    pHeader['FREQUENC']     = ('%d') % (srhRawFits.freqList[frequencyIndex]/1e3 + 0.5)
    pHeader['CDELT1']       = resultArcsecPerPixel
    pHeader['CDELT2']       = resultArcsecPerPixel
    pHeader['CRPIX1']       = srh_x_size // 2
    pHeader['CRPIX2']       = srh_x_size // 2
    pHeader['CTYPE1']       = 'HPLN-TAN'
    pHeader['CTYPE2']       = 'HPLT-TAN'
    pHeader['CUNIT1']       = 'arcsec'
    pHeader['CUNIT2']       = 'arcsec'
    pHeader['PSF_ELLA']     = a # PSF ellipse A arcsec
    pHeader['PSF_ELLB']     = b # PSF ellipse B arcsec
    pHeader['PSF_ELLT']     = ang # PSF ellipse theta deg


    saveFitsIhdu = fits.PrimaryHDU(header=pHeader, data=iImage.astype('float32'))
    saveFitsIpath = saveCleanFitsPath + 'srh_I_%sT%s_%04d.fit'%(srhRawFits.hduList[0].header['DATE-OBS'], fitsTime, srhRawFits.freqList[frequencyIndex]*1e-3 + .5)
    ewLcpPhaseColumn = fits.Column(name='ewLcpPhase', format='D', array = srhRawFits.ewAntPhaLcp[frequencyIndex,:] + srhRawFits.ewLcpPhaseCorrection[frequencyIndex,:])
    ewRcpPhaseColumn = fits.Column(name='ewRcpPhase', format='D', array = srhRawFits.ewAntPhaRcp[frequencyIndex,:] + srhRawFits.ewRcpPhaseCorrection[frequencyIndex,:])
    sLcpPhaseColumn = fits.Column(name='sLcpPhase',   format='D', array = srhRawFits.sAntPhaLcp[frequencyIndex,:] + srhRawFits.sLcpPhaseCorrection[frequencyIndex,:])
    sRcpPhaseColumn = fits.Column(name='sRcpPhase',   format='D', array = srhRawFits.sAntPhaRcp[frequencyIndex,:] + srhRawFits.sRcpPhaseCorrection[frequencyIndex,:])
    saveFitsIExtHdu = fits.BinTableHDU.from_columns([ewLcpPhaseColumn, ewRcpPhaseColumn, sLcpPhaseColumn, sRcpPhaseColumn])
    hduList = fits.HDUList([saveFitsIhdu, saveFitsIExtHdu])
    hduList.writeto(saveFitsIpath, overwrite=True)
    
    saveFitsVhdu = fits.PrimaryHDU(header=pHeader, data=vImage.astype('float32'))
    saveFitsVpath = saveCleanFitsPath + 'srh_V_%sT%s_%04d.fit'%(srhRawFits.hduList[0].header['DATE-OBS'], fitsTime, srhRawFits.freqList[frequencyIndex]*1e-3 + .5)
    hduList = fits.HDUList(saveFitsVhdu)
    hduList.writeto(saveFitsVpath, overwrite=True)

date = '20220314'
calibrate = False
frequency = 8
step = 5
scan = 0
imSize = 1024
cell = 2.45
gainTable = '//home/sergey_lesovoi/SRH0612/20220314/9000_work_1.json'


if (not os.path.exists('cleanMaps/')):
    os.system('mkdir ' + 'cleanMaps/')
if (not os.path.exists('png/')):
    os.system('mkdir ' + 'png/')
if (not os.path.exists('cleanMaps/' + date)):
    os.system('mkdir ' + 'cleanMaps/' + date)
if (not os.path.exists('png/' + date)):
    os.system('mkdir ' + 'png/' + date)
print('Folders create')
    
disk_mask = createDiskMask(imSize, cell, 1100).astype(bool)
disk_model = createDiskMask(imSize, cell, 995)
dL = 2*( 10//2) + 1
kern = NP.ones((dL,dL))
disk_model = scipy.signal.fftconvolve(disk_model,kern) / dL**2
disk_model = disk_model[dL//2:dL//2+imSize,dL//2:dL//2+imSize]
disk_model[disk_model<1e-10] = 0

fitPath = '/home/sergey_lesovoi/SRH0612/20220314'
print(fitPath)
fitNames =  findFits(fitPath,'srh_0612_%s*.fit'%date)
print(fitNames)
fitNames.sort()


ia = IA()
PL.ioff()
PL.figure(figsize = (10,10))

for fileName in fitNames[:]:
    file = SrhFitsFile(fileName, 1024)
    file.useNonlinearApproach = True
    
    for file_scan in range(0,20,step):
        file.getHourAngle(file_scan)
        file.setCalibIndex(file_scan)
        pAngle = file.pAngle
        
        if calibrate:
            file.calibrate(frequency)
            file.vis2uv(file_scan)
            file.centerDisk()
        else:
            file.loadGains(gainTable)

        freq_mhz = int(file.freqList[frequency]/1e3)
        visName = 'MS/' + fileName.split('/')[-1].split('.')[0] + '_' + str(freq_mhz) + '_%03d.ms' % scan
        ms2Table = srh612MS2.Srh612Ms2(visName)
        ms2Table.createMS(file, frequencyChannel = [int(frequency)], phaseCorrect = True, amplitudeCorrect = True)


        imageNameDirty = 'images/%03d_dirty_%d' % (scan, datetime.datetime.now().timestamp())
        imageName = 'images/%03d_%d' % (scan, datetime.datetime.now().timestamp())
                
        tclean(vis = visName,
                imagename = imageNameDirty,
                scan = str(file_scan) + '~' + str(file_scan+4),
                stokes = 'RRLL',
                cell = cell,
                imsize = imSize,
                niter = 0)
        
        casa_im = image(imageNameDirty + '.image')
        im = casa_im.getdata()
        mask = NP.zeros_like(im)
        casa_im.unlock()
        
        mask[0,0] = (im[0,0] > 5000) & disk_mask
        mask[0,1] = (im[0,1] > 5000) & disk_mask
        
        # unmasking eruption
        # mask[0,0,1000:1300,400:600] = 1
        # mask[0,1,1000:1300,400:600] = 1
        
        casa_im = image(imageNameDirty + '.pb')
        casa_im.putdata(mask)
        casa_im.saveas(imageNameDirty + '_mask', overwrite=True)
        casa_im.unlock()
        
        casa_model = image(imageNameDirty + '.model')
        im = casa_model.getdata()
        
        diskTb = file.ZirinQSunTb.getTbAtFrequency(file.freqList[frequency]*1e-6)*1e3
        
        im[0,0] = disk_model * diskTb * file.lm_hd_relation[frequency] / file.convolutionNormCoef
        im[0,1] = disk_model * diskTb * file.lm_hd_relation[frequency] / file.convolutionNormCoef
        
        casa_model.putdata(im)
        casa_model.saveas(imageNameDirty + '_model', overwrite=True)
        casa_model.unlock()
        
        ia.open(imageNameDirty + '.image')
        restoring_beam = ia.restoringbeam()['beams']['*0']['*0']
        a,b,ang = restoring_beam['major']['value'],restoring_beam['minor']['value'],restoring_beam['positionangle']['value']
        ia.close()
                
        tclean(vis = visName,
               imagename = imageName,
               scan = str(file_scan) + '~' + str(file_scan+4),
               stokes = 'RRLL',
               cell = cell,
               imsize = imSize,
               niter = 100000,
               gain = 0.01,
               threshold = 20000,
               deconvolver = 'multiscale',
               scales = [0,2,5,7,10,15],
               usemask = 'user',
               pbcor = False,
               restoringbeam = ['%.2farcsec'%(a*0.8), '%.2farcsec'%(b*0.8), '%.2fdeg'%ang],
               mask = imageNameDirty + '_mask',
               startmodel = imageNameDirty + '_model')
                
        ia.open(imageName + '.image')
        restoring_beam = ia.restoringbeam()#['beams']['*0']['*0']
        rcp = ia.getchunk()[:,:,0,0].transpose()
        lcp = ia.getchunk()[:,:,1,0].transpose()
        ia.close()
        
        ia.open(imageName + '.psf')
        psf = ia.getchunk()[:,:,0,0].transpose()
        ia.close()
        
        # rmtables(tablenames = 'images/%03d*' % scan)
        # rmtables(tablenames = visName)
        
        srh_y_size = imSize
        srh_x_size = imSize
        O = srh_x_size//2
        Q = srh_x_size//4
        scale = AffineTransform(scale=(0.5,0.5))
        shift = AffineTransform(translation=(-srh_y_size/2,-srh_y_size/2))
        rotate = AffineTransform(rotation = -pAngle)
        back_shift = AffineTransform(translation=(srh_y_size/2,srh_y_size/2))
        
        rcp = warp(rcp,(shift + (rotate + back_shift)).inverse)
        lcp = warp(lcp,(shift + (rotate + back_shift)).inverse)
        psf = warp(psf,(shift + (rotate + back_shift)).inverse)
        rcp = warp(rcp,(shift + (scale + back_shift)).inverse)[O-Q:O+Q,O-Q:O+Q]
        lcp = warp(lcp,(shift + (scale + back_shift)).inverse)[O-Q:O+Q,O-Q:O+Q]
            
        saveFitsImages((rcp + lcp)/2, (rcp - lcp)/2, restoring_beam, 'cleanMaps/' + date + '/', file, frequency, file_scan)
 
        PL.clf()
        PL.imshow((rcp + lcp)/2, origin = 'lower', cmap = 'hot', vmin = -5000, vmax = 50000)
        fitsTime = ihhmm_format(file.freqTime[frequency, file_scan])
        PL.savefig('png/%s/i%s_%d_%03d.png'% (date,fitsTime,file.freqList[frequency],scan))

        scan += step
PL.ion()
print('Finish')