#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 09:23:28 2022

@author: sergey_lesovoi
"""
import telebot
from datetime import datetime
import time
#import srhFitsFile612
from py_srh0612_data import Srh0612UVData
from py_srh0306_data import Srh0306UVData
from casatasks import tclean, importuvfits, rmtables
from casatools import image as IA
import pylab as PL
import numpy as NP
from skimage.transform import warp, AffineTransform
from skimage import measure
from srhFitsFile612 import SrhFitsFile
from astropy.io import fits
import ftplib;
import os.path
import os, fnmatch
import glob, shutil
#import argparse

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

def saveFitsImages(iImage, vImage, psfImage, saveCleanFitsPath, srhRawFits, frequencyIndex, array):
    resultArcsecPerPixel = 4.9
    srh_x_size = 1024
    
    contours = (measure.find_contours(psfImage, 0.5))[0]
    con = NP.zeros_like(contours)
    con[:,1] = contours[:,0]
    con[:,0] = contours[:,1]
    sunEll = measure.EllipseModel()
    ellipseEstimateResult = sunEll.estimate(con)
    
    
    pHeader = fits.Header();
    pHeader['DATE-OBS']     = srhRawFits.hduList[0].header['DATE-OBS']
    pHeader['T-OBS']        = srhRawFits.hduList[0].header['TIME-OBS']
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
    if (ellipseEstimateResult):
        pHeader['PSF_ELLA']     = NP.ceil(sunEll.params[2] * resultArcsecPerPixel*10)/10 # PSF ellipse A arcsec
        pHeader['PSF_ELLB']     = NP.ceil(sunEll.params[3] * resultArcsecPerPixel*10)/10 # PSF ellipse B arcsec
        pHeader['PSF_ELLT']     = sunEll.params[4] # PSF ellipse theta rad

    saveFitsIhdu = fits.PrimaryHDU(header=pHeader, data=iImage.astype('float32'))
    saveFitsIpath = saveCleanFitsPath + 'srh_I_%sT%s_%04d.fit'%(srhRawFits.hduList[0].header['DATE-OBS'], srhRawFits.hduList[0].header['TIME-OBS'], srhRawFits.freqList[frequencyIndex]*1e-3 + .5)
    ewLcpPhaseColumn = fits.Column(name='ewLcpPhase', format='D', array = srhRawFits.ewAntPhaLcp[frequencyIndex,:] + srhRawFits.ewLcpPhaseCorrection[frequencyIndex,:])
    ewRcpPhaseColumn = fits.Column(name='ewRcpPhase', format='D', array = srhRawFits.ewAntPhaRcp[frequencyIndex,:] + srhRawFits.ewRcpPhaseCorrection[frequencyIndex,:])
    if (array == '0612'):
        sLcpPhaseColumn = fits.Column(name='sLcpPhase',   format='D', array = srhRawFits.sAntPhaLcp[frequencyIndex,:] + srhRawFits.sLcpPhaseCorrection[frequencyIndex,:])
        sRcpPhaseColumn = fits.Column(name='sRcpPhase',   format='D', array = srhRawFits.sAntPhaRcp[frequencyIndex,:] + srhRawFits.sRcpPhaseCorrection[frequencyIndex,:])
        saveFitsIExtHdu = fits.BinTableHDU.from_columns([ewLcpPhaseColumn, ewRcpPhaseColumn, sLcpPhaseColumn, sRcpPhaseColumn])
    else:
        nLcpPhaseColumn = fits.Column(name='nLcpPhase',   format='D', array = srhRawFits.nAntPhaLcp[frequencyIndex,:] + srhRawFits.nLcpPhaseCorrection[frequencyIndex,:])
        nRcpPhaseColumn = fits.Column(name='nRcpPhase',   format='D', array = srhRawFits.nAntPhaRcp[frequencyIndex,:] + srhRawFits.nRcpPhaseCorrection[frequencyIndex,:])
        saveFitsIExtHdu = fits.BinTableHDU.from_columns([ewLcpPhaseColumn, ewRcpPhaseColumn, nLcpPhaseColumn, nRcpPhaseColumn])
    hduList = fits.HDUList([saveFitsIhdu, saveFitsIExtHdu])
    hduList.writeto(saveFitsIpath)
    
    saveFitsVhdu = fits.PrimaryHDU(header=pHeader, data=vImage.astype('float32'))
    saveFitsVpath = saveCleanFitsPath + 'srh_V_%sT%s_%04d.fit'%(srhRawFits.hduList[0].header['DATE-OBS'], srhRawFits.hduList[0].header['TIME-OBS'], srhRawFits.freqList[frequencyIndex]*1e-3 + .5)
    hduList = fits.HDUList(saveFitsVhdu)
    hduList.writeto(saveFitsVpath)
    return (saveFitsIpath, saveFitsVpath)

bot = telebot.TeleBot('1807793449:AAHtCSQhj0lcAD9jbMUb7R3m8nm5er4bJDw')
srhDate = datetime.now().strftime('%Y%m%d')
srhTime = datetime.now().strftime('%H%M%S')
srhChannel = 0
srhScan = 0
srhAverage = 0
srhArray = '0612'

@bot.message_handler(commands=['start'])
def start_message(message):
    bot.send_message(message.chat.id, 'Hi!\nParameters\ndate ' + srhDate + '\ntime ' + srhTime + '\nchannel ' + 
                     str(srhChannel) + '\nscan ' + str(srhScan) + '\naverage ' + str(srhAverage) + '\narray ' + srhArray)
  
@bot.message_handler(commands=['help'])
def help_message(message):
  bot.send_message(message.chat.id,'Use /prepare_image YYYYMMDD hhmmss f s a A,\nwhere f - frequency channel number, s - scan number, a - averaging, A - array\n' + 
                   'or /setDate YYYYMMDD\n' + 
                   '/setTime hhmmss\n' + 
                   '/setChannel f\n' + 
                   '/setScan s\n' + 
                   '/setAverage a\n' + 
                   '/setArray LLHH\n' + 
                   'For example: prepare_image 20220113 031010 2 0 3 0612\n')
  
# def stop(message):
#   if message.from_user.username == cfg.Father:	
#     pid = str(os.getpid())
#     stoper = open('ozerx/stoper.bat', 'w')
#     stoper.write("Taskkill /PID " + pid + " /F")
#     stoper.close()
#     os.system('C:/Users/smp/Desktop/SMP/ozerx/stoper.bat')
#   else:
#     bot.send_message(message.chat.id, "Ты не Создатель бота; у тебя нет админ-прав, проваливай!

@bot.message_handler(commands=['setDate'])
def setDate_message(message):
    global srhDate
    args = message.text.split()[1:]
    if (len(args) > 0):
        srhDate = args[0]
    bot.send_message(message.chat.id, 'Parameters\ndate ' + srhDate + '\ntime ' + srhTime + '\nchannel ' + 
                     str(srhChannel) + '\nscan ' + str(srhScan) + '\naverage ' + str(srhAverage) + '\narray ' + srhArray)

@bot.message_handler(commands=['setTime'])
def setTime_message(message):
    global srhTime
    args = message.text.split()[1:]
    if (len(args) > 0):
        srhTime = args[0]
    bot.send_message(message.chat.id, 'Parameters\ndate ' + srhDate + '\ntime ' + srhTime + '\nchannel ' + 
                     str(srhChannel) + '\nscan ' + str(srhScan) + '\naverage ' + str(srhAverage) + '\narray ' + srhArray)
                     
@bot.message_handler(commands=['setChannel'])
def setChannel_message(message):
    global srhChannel
    args = message.text.split()[1:]
    if (len(args) > 0):
        srhChannel = args[0]
    bot.send_message(message.chat.id, 'Parameters\ndate ' + srhDate + '\ntime ' + srhTime + '\nchannel ' + 
                     str(srhChannel) + '\nscan ' + str(srhScan) + '\naverage ' + str(srhAverage) + '\narray ' + srhArray)
                     
@bot.message_handler(commands=['setScan'])
def setScan_message(message):
    global srhScan
    args = message.text.split()[1:]
    if (len(args) > 0):
        srhScan = args[0]
    bot.send_message(message.chat.id, 'Parameters\ndate ' + srhDate + '\ntime ' + srhTime + '\nchannel ' + 
                     str(srhChannel) + '\nscan ' + str(srhScan) + '\naverage ' + str(srhAverage) + '\narray ' + srhArray)

@bot.message_handler(commands=['setAverage'])
def setAverage_message(message):
    global srhAverage
    args = message.text.split()[1:]
    if (len(args) > 0):
        srhAverage = args[0]
    bot.send_message(message.chat.id, 'Parameters\ndate ' + srhDate + '\ntime ' + srhTime + '\nchannel ' + 
                     str(srhChannel) + '\nscan ' + str(srhScan) + '\naverage ' + str(srhAverage) + '\narray ' + srhArray)
                     
@bot.message_handler(commands=['setArray'])
def setArray_message(message):
    global srhArray
    args = message.text.split()[1:]
    if (len(args) > 0):
        srhArray = args[0]
    bot.send_message(message.chat.id, 'Parameters\ndate ' + srhDate + '\ntime ' + srhTime + '\nchannel ' + 
                     str(srhChannel) + '\nscan ' + str(srhScan) + '\naverage ' + str(srhAverage) + '\narray ' + srhArray)
                     
@bot.message_handler(commands=['prepare_image'])
def prepare_image_message(message):
    args = message.text.split()
    if (len(args) > 1):
        imageDate = args[1]
    else:
        imageDate = srhDate
    if (len(args) > 2):
        imageTime = args[2]
    else:
        imageTime = srhTime
    if (len(args) > 3):
        imageChannel = int(args[3])
    else:
        imageChannel = int(srhChannel)
    if (len(args) > 4):
        imageScan = int(args[4])
    else:
        imageScan = int(srhScan)
    if (len(args) > 5):
        imageAverage = int(args[5])
    else:
        imageAverage = int(srhAverage)
    if (len(args) > 6):
        imageArray = args[6]
    else:
        imageArray = srhArray
        
    if (imageArray == '0612'):
        srhImageName = 'srh_' + imageArray + '_' + imageDate + 'T' + imageTime
    else:
        srhImageName = 'srh_' + imageDate + 'T' + imageTime
    srhFitName = srhImageName + '.fit'
    srhUvFitsName = srhImageName + '.uvfits'
    srhVisName = srhImageName + '.ms'
    srhFitPath = '../SRH_DATA/SRH/SRH' + imageArray + '/' + imageDate + '/'
    srhFitPathName = srhFitPath + srhFitName
    if (not os.path.exists(srhFitPathName)):
        bot.send_message(message.chat.id, srhFitPathName + ' not found')
        if (imageArray == '0612'):
            existingFitNames = findFits(srhFitPath,'srh_' + imageArray + '_' + imageDate + 'T' + imageTime[0:3] + '*.fit')
        else:
            existingFitNames = findFits(srhFitPath,'srh_'                    + imageDate + 'T' + imageTime[0:3] + '*.fit')
#            existingFitNames = findFits('../SRH_DATA/SRH/SRH0306/20220302/','*.fit')
        if (len(existingFitNames) > 0):
            existingFitNames.sort()
            srhFitPathName = existingFitNames[0]
            bot.send_message(message.chat.id, 'Trying ' + srhFitPathName)
        else:
            bot.send_message(message.chat.id, 'findFits failure ' + srhFitPath + 'srh_'                    + imageDate + 'T' + imageTime[0:4] + '*.fit')

    if (imageArray == '0306'):
        srhUV = Srh0306UVData()
    else:
        srhUV = Srh0612UVData()
        
    try:
        srhUV.write_uvfits(srhFitPathName,srhUvFitsName,calibrating=True,frequency=imageChannel,scan=imageScan, averaging=imageAverage)
        bot.send_message(message.chat.id, 'fits:' + srhFitPathName + ' --> uvfits:' + srhUvFitsName)
        rmtables(srhVisName)
        importuvfits(fitsfile = srhUvFitsName, vis = 'botMS/' + srhVisName)
        tclean(vis='botMS/' + srhVisName, imagename='botImages/' + srhImageName,imsize=1024, cell=['3arcsec'],niter=10000, stokes = 'RRLL', threshold = '30000Jy')
        ia = IA()
        ia.open('botImages/' + srhImageName + '.image')
        rcp = ia.getchunk()[:,:,0,0].transpose()
        lcp = ia.getchunk()[:,:,1,0].transpose()
        ia.close()
        
        ia.open('botImages/' + srhImageName + '.psf')
        psf = ia.getchunk()[:,:,0,0].transpose()
        ia.close()
        
        file = SrhFitsFile(srhFitPathName, 1024)
        file.useNonlinearApproach = True
        file.getHourAngle(0)
        pAngle = file.pAngle

        srh_y_size = 1024
        srh_x_size = 1024
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
        psf = warp(psf,(shift + (scale + back_shift)).inverse)[O-Q:O+Q,O-Q:O+Q]
            
        PL.imshow(rcp+lcp,origin='lower',cmap='jet',vmin=0,vmax=1e5)
        PL.title('SRH %s, %s, %.2f GHz'%(imageDate, imageTime, srhUV.freq_array*1e-9))
        PL.axis('off')
        PL.savefig('botImages/tempI.png')
        fImage = open('botImages/tempI.png','rb')
        bot.send_photo(message.chat.id, fImage)
        fImage.close()
        os.remove('botImages/tempI.png')
        
        savedIpath, savedVpath = saveFitsImages((rcp + lcp)/2, (rcp - lcp)/2, psf, 'botCleanMaps/', file, imageChannel, srhArray)
        store = ftplib.FTP('ftp.rao.istp.ac.ru','ftpwriter','cn5HeuA4u')
        fi = open(savedIpath,'rb')
        store.storbinary('STOR /SRH/botCleanMaps/' + savedIpath.split('/')[-1], fi)
        fi.close()
        fi = open(savedVpath,'rb')
        store.storbinary('STOR /SRH/botCleanMaps/' + savedVpath.split('/')[-1], fi)
        fi.close()
        store.close()
        bot.send_message(message.chat.id,'I data were saved at ftp://ftp.rao.istp.ac.ru/SRH/botCleanMaps/' + savedIpath.split('/')[-1])
        bot.send_message(message.chat.id,'V data were saved at ftp://ftp.rao.istp.ac.ru/SRH/botCleanMaps/' + savedVpath.split('/')[-1])
        
        botMS = glob.glob("botMS/*")
        for ms in botMS:
            shutil.rmtree(ms)
        botImages = glob.glob("botImages/*")
        for images in botImages:
            shutil.rmtree(images)
        botCleanMaps = glob.glob("botCleanMaps/*")
        for cleanMap in botCleanMaps:
            os.remove(cleanMap)
    except IOError:
        bot.send_message(message.chat.id, 'Sorry, there is none any image at this time ' + imageTime)
        
        
@bot.message_handler(content_types=["text"])
def repeat_all_messages(message): 
    if (message.text[0] == '/'):
    	bot.send_message(message.chat.id, 'Unknown command')
    else:
    	bot.send_message(message.chat.id, message.text)
#------------------------------------------------------------------------------  
while True:
	try:
		bot.polling(none_stop = True)
	except Exception as msg:
		print(msg)
		time.sleep(5)

