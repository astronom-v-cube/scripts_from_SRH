# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 01:47:01 2016

@author: Sergey
"""
import os, fnmatch;
from srhFitsFile36 import SrhFitsFile
import numpy;
import datetime as DT;
import pylab as LAB;
from PyQt5 import QtWidgets,QtCore;
import sys;
import ftplib;
from optparse import OptionParser;
import matplotlib;
from scipy.io import readsav

from astropy.io import fits

dt_major = 3600.;
dt_minor = 900.;
freq_first = 2;
freq_intervals = 5;

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

def findFits(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

application = QtWidgets.QApplication.instance()
if not application:    
    application = QtWidgets.QApplication(sys.argv)


#os.environ['DISPLAY'] = ':0.0';

dateName = DT.datetime.now().strftime("%Y%m%d");
#fitNames = findFits('/home/svlesovoi/SRH_DATA/SRH36/' + dateName ,'*.fit')
fitNames = findFits('/home/svlesovoi/Documents/Python Scripts/SRH36_temp_20210322_2','*.fit')
    
fitNames.sort();
fitNames = fitNames[1:];
fdSrhImageName = open('/home/svlesovoi/SRH/srh_currentFileName.txt','w');
fdSrhImageName.write(fitNames[-1]);
fdSrhImageName.close();
visFirst = 32*0
visLast = 32*16 - 1

for fName in fitNames:
    sF = SrhFitsFile(fName, 256)
    if (fName == fitNames[0]):
        freqAmount  = sF.freqListLength;
        freqList = sF.freqList;
        startTime= DT.datetime.strptime(sF.hduList[0].header['DATE-OBS'] + ' ' + sF.hduList[0].header['TIME-OBS'], '%Y-%m-%d %H:%M:%S');
        times = sF.freqTime;
        r_real = abs(sF.visRcp.real[:,:,0:3007])# / float(sF.hduList[0].header['VIS_MAX'])
        r_imag = abs(sF.visRcp.imag[:,:,0:3007])#/ float(sF.hduList[0].header['VIS_MAX'])
        fluxRcp = numpy.mean(numpy.sqrt(r_real**2. + r_imag**2.),axis=2);
        r_real = abs(sF.visLcp.real[:,:,0:3007])#/ float(sF.hduList[0].header['VIS_MAX'])
        r_imag = abs(sF.visLcp.imag[:,:,0:3007])#/ float(sF.hduList[0].header['VIS_MAX'])
        fluxLcp = numpy.mean(numpy.sqrt(r_real**2. + r_imag**2.),axis=2);
        r_real = 0;
        r_imag = 0;
        
        
    else:
        if sF.freqTime[0,1] - times[0,-1] > 600.:
            insTime = numpy.zeros((freqAmount,2),dtype=float);
            insRcp = numpy.zeros((freqAmount,2),dtype=float);
            insLcp = numpy.zeros((freqAmount,2),dtype=float);
            for fff in range(freqAmount):
                insTime[fff,0] = times[fff,-1];
                insTime[fff,1] = sF.freqTime[fff,0]
            times = numpy.concatenate((times,insTime),axis=1);
            fluxRcp = numpy.concatenate((fluxRcp,insRcp),axis=1);
            fluxLcp = numpy.concatenate((fluxLcp,insLcp),axis=1);
        times = numpy.concatenate((times,sF.freqTime),axis=1);
        r_real = abs(sF.visRcp.real[:,:,0:3007])#/ float(sF.hduList[0].header['VIS_MAX'])
        r_imag = abs(sF.visRcp.imag[:,:,0:3007])#/ float(sF.hduList[0].header['VIS_MAX'])
        fluxRcp = numpy.concatenate((fluxRcp, numpy.mean(numpy.sqrt(r_real**2. + r_imag**2.),axis=2)),axis=1);
        r_real = abs(sF.visLcp.real[:,:,0:3007])#/ float(sF.hduList[0].header['VIS_MAX'])
        r_imag = abs(sF.visLcp.imag[:,:,0:3007])#/ float(sF.hduList[0].header['VIS_MAX'])
        fluxLcp = numpy.concatenate((fluxLcp, numpy.mean(numpy.sqrt(r_real**2. + r_imag**2.),axis=2)),axis=1);
        r_real = 0;
        r_imag = 0;
        
    sF.close();

vNorm = numpy.mean(fluxRcp[:,:],axis=1) / numpy.mean(fluxLcp[:,:],axis=1);

c_list = matplotlib.colors.LinearSegmentedColormap.from_list(LAB.cm.datad['gist_rainbow'], colors=['r','g','b'], N = freqAmount);
saveFreqs = numpy.linspace(0,freqAmount-1,freqAmount,dtype=numpy.uint8)


# creating corrPlot

#LAB.clf();
fig = LAB.figure(figsize = (16,8));
fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
fig.tight_layout()
commonTitle = 'Correlation plot ' + startTime.strftime('%Y %B %d') + '. Stokes I, V at 3-6 GHz';
fig.suptitle(commonTitle,fontsize=8);

sub = fig.add_subplot(1,1,1);
sub.set_ylabel('correlation coefficient');
sub.set_xlabel('UT');
sub.xaxis.set_major_locator(LAB.MultipleLocator(dt_major));
sub.xaxis.set_major_formatter(LAB.FuncFormatter(hhmm_format));
sub.xaxis.set_minor_locator(LAB.MultipleLocator(dt_minor));
sub.set_xlim((0.*3600.,11.*3600.));
sub.set_ylim((-0.002,0.01));

saveTime = numpy.zeros((saveFreqs.shape[0],times.shape[1]));
saveI = numpy.zeros((saveFreqs.shape[0],times.shape[1]));
saveV = numpy.zeros((saveFreqs.shape[0],times.shape[1]));
for ff in saveFreqs:
    sub.scatter(times[ff,:],0.5*(fluxRcp[ff,:] + vNorm[ff]*fluxLcp[ff,:]),color=c_list(ff), s=0.5, linewidths = 0)
    sub.scatter(times[ff,:],0.5*(fluxRcp[ff,:] - vNorm[ff]*fluxLcp[ff,:]),color=c_list(ff), s=0.5, linewidths = 0)
    saveTime[ff,:] = times[ff,:];
    saveI[ff,:] = 0.5*(fluxRcp[ff,:] + vNorm[ff]*fluxLcp[ff,:]);
    saveV[ff,:] = 0.5*(fluxRcp[ff,:] - vNorm[ff]*fluxLcp[ff,:]);
corrPlotName = 'fCorrPlot'+ dateName;
corrPlotName_png = corrPlotName + '.png';
#LAB.savefig(corrPlotName_png);        
    

# writing srh_cp file
    
#dataFormat = str(times.shape[1]) + 'D';
#freqsColumn = fits.Column(name='frequencies',format='D',array=freqList[saveFreqs]);
#timeColumn = fits.Column(name='time',format=dataFormat,array=saveTime);
#IColumn = fits.Column(name='I',format=dataFormat,array=saveI);
#VColumn = fits.Column(name='V',format=dataFormat,array=saveV);
#
#fTableHdu = fits.BinTableHDU.from_columns([freqsColumn]);
#dTableHdu = fits.BinTableHDU.from_columns([timeColumn, IColumn, VColumn]);
#pHeader = fits.Header();
#pHeader['DATE-OBS']     = sF.hduList[0].header['DATE-OBS'];
#pHeader['TIME-OBS']     = sF.hduList[0].header['TIME-OBS'];
#pHeader['INSTRUME']     = 'SRH'
#pHeader['ORIGIN']       = 'ISTP'
#pHeader['OBS-LAT']      = '51.759'
#pHeader['OBS-LONG']     = '102.217'
#pHeader['OBS-ALT']      = '799'
#pHeader['FR_CHAN']      = '10'
#pHdu = fits.PrimaryHDU(header=pHeader);
#hduList = fits.HDUList([pHdu, fTableHdu, dTableHdu]);
#hduList.writeto('srh_cp_' + dateName + '.fits',clobber=True);
#hduList.close();
#


# sending files via FTP

#fd = ftplib.FTP('192.168.0.1','sergey','jill21');
#
#fdCPlotName = open('corrPlotName.txt','w');
#fdCPlotName.write(corrPlotName_png);
#fdCPlotName.close();
#fi = open(corrPlotName_png,'rb');
#fd.storbinary('STOR /Public/sergey/corrPlots/' + corrPlotName_png,fi);
#fi.close();
#fi = open('corrPlotName.txt','rb');
#fd.storbinary('STOR /Public/sergey/corrPlots/corrPlotName.txt',fi);
#fi.close();
#fi = open('srh_cp_' + dateName + '.fits','rb');
#fd.storbinary('STOR /Public/sergey/corrPlots/' + 'srh_cp_' + dateName + '.fits',fi);
#fi.close();
#
#srhImageNameIDL = 'srh_'+ dateName + '.png';
#srhImageName_png = 'srh_'+ DT.datetime.now().strftime("%Y%m%d_%H%M") + '.png';
#os.rename(srhImageNameIDL,srhImageName_png)
#fi = open(srhImageName_png,'rb');
#fd.storbinary('STOR /Public/sergey/corrPlots/' + srhImageName_png,fi);
#fi.close();
#
#fdSrhImageName = open('srhImageName.txt','w');
#fdSrhImageName.write(srhImageName_png);
#fdSrhImageName.close();
#fi = open('srhImageName.txt','rb');
#fd.storbinary('STOR /Public/sergey/corrPlots/srhImageName.txt',fi);
#fi.close();
#fd.quit();
