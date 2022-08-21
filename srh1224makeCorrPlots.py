#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 06:32:14 2021

@author: svlesovoi
"""

import datetime as DT;
from astropy.io import fits
import srh0612Utils
import numpy as NP
import pylab as PL
from srhFitsFile1224 import SrhFitsFile
import matplotlib
import ftplib;
from ZirinTb import ZirinTb

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

dateName = DT.datetime.now().strftime("%Y%m%d")
try:
    previousCpFits = fits.open('srh_1224_cp_' + dateName + '.fits')
except FileNotFoundError:
    previousCpFits = 0

dt_major = 3600.;
dt_minor = 900.;

currentRawFitses = srh0612Utils.findFits('SRH1224/' + dateName ,'*.fit')  
#currentRawFitses = srh0612Utils.findFits('../SRH_1224/' + dateName ,'*.fit')  
currentRawFitses.sort()
currentRawFitses = NP.array(currentRawFitses[1:])

Zi = ZirinTb()

if (previousCpFits):
    previousRawFitses = previousCpFits[3].data['rawFitsNames']
    fitNames = currentRawFitses[previousRawFitses.shape[0]:]
else:
    fitNames = currentRawFitses

for fName in fitNames:
    print(fName)
    sF = SrhFitsFile(fName, 256)
    if (fName == fitNames[0]):
        freqAmount  = sF.freqListLength;
        freqList = sF.freqList;
        startTime= DT.datetime.strptime(sF.hduList[0].header['DATE-OBS'] + ' ' + sF.hduList[0].header['TIME-OBS'], '%Y-%m-%d %H:%M:%S');
        times = sF.freqTime;
        r_real = abs(sF.visRcp.real[:,:,0:8840])
        r_imag = abs(sF.visRcp.imag[:,:,0:8840])
        corrFluxRcp = NP.mean(NP.sqrt(r_real**2. + r_imag**2.),axis=2);
        r_real = abs(sF.visLcp.real[:,:,0:8840])
        r_imag = abs(sF.visLcp.imag[:,:,0:8840])
        corrFluxLcp = NP.mean(NP.sqrt(r_real**2. + r_imag**2.),axis=2);
        r_real = 0;
        r_imag = 0;
        ampFluxRcp = NP.mean(sF.ampRcp[:,:,:], axis = 2)
        ampFluxLcp = NP.mean(sF.ampLcp[:,:,:], axis = 2)
    else:
        if sF.freqTime[0,1] - times[0,-1] > 600.:
            insTime = NP.zeros((freqAmount,2),dtype=float);
            insRcp = NP.zeros((freqAmount,2),dtype=float);
            insLcp = NP.zeros((freqAmount,2),dtype=float);
            for fff in range(freqAmount):
                insTime[fff,0] = times[fff,-1];
                insTime[fff,1] = sF.freqTime[fff,0]
            times = NP.concatenate((times,insTime),axis=1);
            corrFluxRcp = NP.concatenate((corrFluxRcp,insRcp),axis=1);
            corrFluxLcp = NP.concatenate((corrFluxLcp,insLcp),axis=1);
        times = NP.concatenate((times,sF.freqTime),axis=1);
        r_real = abs(sF.visRcp.real[:,:,0:8840])
        r_imag = abs(sF.visRcp.imag[:,:,0:8840])
        corrFluxRcp = NP.concatenate((corrFluxRcp, NP.mean(NP.sqrt(r_real**2. + r_imag**2.),axis=2)),axis=1);
        r_real = abs(sF.visLcp.real[:,:,0:8840])
        r_imag = abs(sF.visLcp.imag[:,:,0:8840])
        corrFluxLcp = NP.concatenate((corrFluxLcp, NP.mean(NP.sqrt(r_real**2. + r_imag**2.),axis=2)),axis=1);
        r_real = 0;
        r_imag = 0;
        ampFluxRcp = NP.concatenate((ampFluxRcp, NP.mean(sF.ampRcp[:,:,:], axis = 2)), axis = 1)
        ampFluxLcp = NP.concatenate((ampFluxLcp, NP.mean(sF.ampLcp[:,:,:], axis = 2)), axis = 1)
        
    sF.close();
    
vNorm = NP.array([1.00220427, 1.00607914, 1.01583696, 1.009751  , 1.00277865,
       0.99725935, 0.99515777, 1.00312649, 1.0008489 , 1.00158222,
       1.00313546, 1.0059643 , 1.01301051, 0.99926368, 0.99368663,
       0.99481965])
    
saveFreqs = NP.linspace(0,freqAmount-1,freqAmount,dtype=NP.uint8)
saveTime = NP.zeros((saveFreqs.shape[0],times.shape[1]))
saveCorrI = NP.zeros((saveFreqs.shape[0],times.shape[1]))
saveCorrV = NP.zeros((saveFreqs.shape[0],times.shape[1]))
saveFluxI = NP.zeros((saveFreqs.shape[0],ampFluxRcp.shape[1]))
saveFluxV = NP.zeros((saveFreqs.shape[0],ampFluxRcp.shape[1]))

for ff in saveFreqs:
    ampFluxRcp[ff,:] *= vNorm[ff]

for ff in saveFreqs:
    saveTime[ff,:] = times[ff,:]
    saveCorrI[ff,:] = 0.5*(corrFluxRcp[ff,:] + corrFluxLcp[ff,:])
    saveCorrV[ff,:] = 0.5*(corrFluxRcp[ff,:] - corrFluxLcp[ff,:])
    saveFluxI[ff,:] = 0.5*(ampFluxRcp[ff,:] + ampFluxLcp[ff,:])
    saveFluxV[ff,:] = 0.5*(ampFluxRcp[ff,:] - ampFluxLcp[ff,:])

try:
    previousCpZerosFits = fits.open('srh_1224_cp_zeros.fits')
    corrZeros = previousCpZerosFits[2].data['corrI']
    fluxZeros = previousCpZerosFits[2].data['fluxI']
except FileNotFoundError:
    corrZeros = NP.zeros_like(freqList)
    fluxZeros = NP.zeros_like(freqList)

try:
    previousCpFluxNormFits = fits.open('srh_1224_cp_fluxNorm.fits')
    fluxNormI = previousCpFluxNormFits[2].data['fluxNormI']
except FileNotFoundError:
    fluxNormI = NP.ones_like(freqList)

for ff in saveFreqs:
#    saveCorrI[ff,:] -= corrZeros[ff]
    saveFluxI[ff,:] -= fluxZeros[ff]
    saveFluxI[ff,:] *= fluxNormI[ff]
    saveFluxV[ff,:] *= fluxNormI[ff]

if (previousCpFits):
    previousTime = previousCpFits[2].data['time']
    previousCorrI = previousCpFits[2].data['I']
    previousCorrV = previousCpFits[2].data['V']
    previousFluxI = previousCpFits[2].data['flux_I']
    previousFluxV = previousCpFits[2].data['flux_V']
    saveTime = NP.concatenate((previousTime, saveTime), axis=1)
    saveCorrI = NP.concatenate((previousCorrI, saveCorrI),axis=1)
    saveCorrV = NP.concatenate((previousCorrV, saveCorrV),axis=1)
    saveFluxI = NP.concatenate((previousFluxI, saveFluxI),axis=1)
    saveFluxV = NP.concatenate((previousFluxV, saveFluxV),axis=1)
    
# creating corrPlot
c_list = matplotlib.colors.LinearSegmentedColormap.from_list(PL.cm.datad['gist_rainbow'], colors=['r','g','b'], N = freqAmount)
fig = PL.figure(figsize = (16,8));
fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.05)
fig.tight_layout()
commonTitle = 'Correlation plot ' + startTime.strftime('%Y %B %d') + '. Stokes I, V at 12-24 GHz'
fig.suptitle(commonTitle,fontsize=14)

sub = fig.add_subplot(1,1,1)
sub.set_ylabel('correlation coefficient')
sub.set_xlabel('UT')
sub.xaxis.set_major_locator(PL.MultipleLocator(dt_major))
sub.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
sub.xaxis.set_minor_locator(PL.MultipleLocator(dt_minor))
sub.set_xlim((0.*3600.,10.*3600.))
sub.set_ylim((-0.0005,0.005))
for ff in saveFreqs:
    sub.scatter(saveTime[ff,:],saveCorrI[ff,:],color=c_list(ff), s=0.5, linewidths = 0,label='%d MHz'%(freqList[ff]*1e-3))
    sub.scatter(saveTime[ff,:],saveCorrV[ff,:],color=c_list(ff), s=0.5, linewidths = 0)
    sub.legend(markerscale=10)
corrPlotName = 'fCorrPlot_1224_'+ dateName
corrPlotName_png = corrPlotName + '.png'
PL.savefig(corrPlotName_png)

# creating fluxPlot
fig = PL.figure(figsize = (16,8))
fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
fig.tight_layout()
commonTitle = 'Flux plot ' + startTime.strftime('%Y %B %d') + '. Stokes I, V at 12-24 GHz'
fig.suptitle(commonTitle,fontsize=14)

sub = fig.add_subplot(1,1,1)
sub.set_ylabel('s.f.u.')
sub.set_xlabel('UT')
sub.xaxis.set_major_locator(PL.MultipleLocator(dt_major))
sub.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
sub.xaxis.set_minor_locator(PL.MultipleLocator(dt_minor))
sub.set_xlim((0.*3600.,10.*3600.))
sub.set_ylim((-2e1, 5e3))

for ff in saveFreqs:
    sub.scatter(saveTime[ff,:],saveFluxI[ff,:],color=c_list(ff), s=1.5, linewidths = 0,label='%d MHz'%(freqList[ff]*1e-3))
    sub.scatter(saveTime[ff,:],saveFluxV[ff,:],color=c_list(ff), s=1.5, linewidths = 0)
    sub.legend(markerscale=5)
fluxPlotName = 'fFluxPlot_1224_'+ dateName
fluxPlotName_png = fluxPlotName + '.png'
PL.savefig(fluxPlotName_png)

# writing srh_cp file
    
dataFormat = str(saveTime.shape[1]) + 'D'
freqsColumn = fits.Column(name='frequencies',format='D',array=freqList[saveFreqs])
timeColumn = fits.Column(name='time',format=dataFormat,array=saveTime)
IColumn = fits.Column(name='I',format=dataFormat,array=saveCorrI)
VColumn = fits.Column(name='V',format=dataFormat,array=saveCorrV)
IAmpColumn = fits.Column(name='flux_I',format=dataFormat,array=saveFluxI)
VAmpColumn = fits.Column(name='flux_V',format=dataFormat,array=saveFluxV)

fTableHdu = fits.BinTableHDU.from_columns([freqsColumn]);
dTableHdu = fits.BinTableHDU.from_columns([timeColumn, IColumn, VColumn, IAmpColumn, VAmpColumn])
pHeader = fits.Header()
pHeader['DATE-OBS']     = sF.hduList[0].header['DATE-OBS']
pHeader['TIME-OBS']     = sF.hduList[0].header['TIME-OBS']
pHeader['INSTRUME']     = 'SRH'
pHeader['ORIGIN']       = 'ISTP'
pHeader['OBS-LAT']      = '51.759'
pHeader['OBS-LONG']     = '102.217'
pHeader['OBS-ALT']      = '799'
pHeader['FR_CHAN']      = '10'

fitsNamesColumn = fits.Column(name='rawFitsNames', format='A256', array=currentRawFitses)
aTableHdu = fits.TableHDU.from_columns([fitsNamesColumn])

pHdu = fits.PrimaryHDU(header=pHeader);
hduList = fits.HDUList([pHdu, fTableHdu, dTableHdu, aTableHdu]);
hduList.writeto('srh_1224_cp_' + dateName + '.fits',clobber=True);
hduList.close();


#zeroCorrI = NP.array(saveCorrI[:,90:120].mean(axis=1))
#zeroFluxI = NP.array(saveFluxI[:,90:120].mean(axis=1))
#
#zeroCorrIColumn = fits.Column(name='corrI',format='D',array=zeroCorrI)
#zeroFluxIColumn = fits.Column(name='fluxI',format='D',array=zeroFluxI)
#
#dTableHdu = fits.BinTableHDU.from_columns([zeroCorrIColumn, zeroFluxIColumn])
#pHdu = fits.PrimaryHDU(header=pHeader);
#hduList = fits.HDUList([pHdu, fTableHdu, dTableHdu]);
#hduList.writeto('srh_1224_cp_zeros_' + dateName + '.fits',clobber=True);
#hduList.close();
#for ff in saveFreqs:
#    saveFluxI[ff] -= zeroFluxI[ff]
#
#fluxINorm = 1/saveFluxI[:,0:50].mean(axis=1) * Zi.getSfuAtFrequency(freqList*1e-6)
#normFluxIColumn = fits.Column(name='fluxNormI',format='D',array=fluxINorm)
#dTableHdu = fits.BinTableHDU.from_columns([normFluxIColumn])
#pHdu = fits.PrimaryHDU(header=pHeader);
#hduList = fits.HDUList([pHdu, fTableHdu, dTableHdu]);
#hduList.writeto('srh_1224_cp_fluxNorm_' + dateName + '.fits',clobber=True);
#hduList.close();
#for ff in saveFreqs:
#    saveFluxI[ff] *= fluxINorm[ff]

fd = ftplib.FTP('10.1.1.9','sergey','jill21');

fdCPlotName = open('corrPlotName.txt','w');
fdCPlotName.write(corrPlotName_png);
fdCPlotName.close();
fi = open(corrPlotName_png,'rb');
fd.storbinary('STOR /Public/sergey/corrPlots/' + corrPlotName_png,fi);
fi.close();
fi = open('corrPlotName.txt','rb');
fd.storbinary('STOR /Public/sergey/corrPlots/corrPlotName.txt',fi);
fi.close();
fi = open('srh_1224_cp_' + dateName + '.fits','rb');
fd.storbinary('STOR /Public/sergey/corrPlots/' + 'srh_1224_cp_' + dateName + '.fits',fi);
fi.close();

fdFPlotName = open('corrMapName.txt','w');
fdFPlotName.write(fluxPlotName_png);
fdFPlotName.close();
fi = open(fluxPlotName_png,'rb');
fd.storbinary('STOR /Public/sergey/corrPlots/' + fluxPlotName_png,fi);
fi.close();
fi = open('corrMapName.txt','rb');
fd.storbinary('STOR /Public/sergey/corrPlots/corrMapName.txt',fi);
fi.close();
