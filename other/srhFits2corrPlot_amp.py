# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 01:47:01 2016

@author: Sergey
"""
import os, fnmatch;
import srhFitsFile;
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

year = 2020
month = 9
day = 13
dateName = DT.datetime(year,month,day).strftime("%Y%m%d");
#fitNames = findFits(('/home/svlesovoi/SRH/%4d/%02d/%02d'%(year,month,day)),'mf_' + dateName + '*.fit')
fitNames = findFits(('/var/run/media/svlesovoi/SRH_DATA/SRH/%4d/%02d/%02d'%(year,month,day)),'mf_' + dateName + '*.fit')
#fitNames = findFits(('/home/sergeylesovoi/SRH/%4d/%02d/%02d'%(year,month,day)),'mf_' + dateName + '*.fit')

fitNames.sort();
fitNames = fitNames[1:];
#fdSrhImageName = open('/home/serg/SRH/mf_currentFileName.txt','w');
#fdSrhImageName.write(fitNames[-1]);
#fdSrhImageName.close();
visFirst = 32*0
visLast = 32*16 - 1

for fName in fitNames:
    sF = srhFitsFile.SrhFitsFile(fName, 256)
    if (fName == fitNames[0]):
        freqAmount  = sF.freqListLength;
        freqList = sF.freqList;
        startTime= DT.datetime.strptime(sF.hduList[0].header['DATE-OBS'] + ' ' + sF.hduList[0].header['TIME-OBS'], '%Y/%m/%d %H:%M:%S.%f');
        times = sF.freqTime;
        r_real = abs(sF.visRcp.real[:,:,0:512])
        r_imag = abs(sF.visRcp.imag[:,:,0:512])
        fluxRcp = numpy.mean(numpy.sqrt(r_real**2. + r_imag**2.),axis=2);
        r_real = abs(sF.visLcp.real[:,:,0:512])
        r_imag = abs(sF.visLcp.imag[:,:,0:512])
        fluxLcp = numpy.mean(numpy.sqrt(r_real**2. + r_imag**2.),axis=2);
        r_real = 0;
        r_imag = 0;
        
        rAmp = sF.ampRcp
        lAmp = sF.ampLcp
        rVis = sF.visRcp[:,:,512:602]
        lVis = sF.visLcp[:,:,512:602]
        antA = sF.antennaA[512:602]
        antB = sF.antennaB[512:602]
        ant = sF.antennaA[800:848]
        
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
        r_real = abs(sF.visRcp.real[:,:,visFirst:visLast])
        r_imag = abs(sF.visRcp.imag[:,:,visFirst:visLast])
        fluxRcp = numpy.concatenate((fluxRcp, numpy.mean(numpy.sqrt(r_real**2. + r_imag**2.),axis=2)),axis=1);
        r_real = abs(sF.visLcp.real[:,:,visFirst:visLast])
        r_imag = abs(sF.visLcp.imag[:,:,visFirst:visLast])
        fluxLcp = numpy.concatenate((fluxLcp, numpy.mean(numpy.sqrt(r_real**2. + r_imag**2.),axis=2)),axis=1);
        r_real = 0;
        r_imag = 0;
        
        rAmp = numpy.concatenate((rAmp, sF.ampRcp), axis = 1)
        lAmp = numpy.concatenate((lAmp, sF.ampLcp), axis = 1)
        rVis = numpy.concatenate((rVis, sF.visRcp[:,:,512:602]), axis = 1)
        lVis = numpy.concatenate((lVis, sF.visLcp[:,:,512:602]), axis = 1)
        
    sF.close();

vNorm = numpy.mean(fluxRcp[:,:],axis=1) / numpy.mean(fluxLcp[:,:],axis=1);

#fdNorpData = open('norp_data.txt','rb');
#norpLines = fdNorpData.readlines();
#norpData = [];
#for i in norpLines:
#	norpData.append([float(j) for j in i.split()]);
#norpData = numpy.array(norpData);
#fdNorpData.close();
#
#norp_freq = norpData[:,0]*1e3;
#norp_flux = norpData[:,1];
#norp_fit_coef = numpy.polyfit(norp_freq,norp_flux,2);
#
#fit_freq = 10e3*numpy.arange(1000.)/1000.;
#fit_flux = numpy.polyval(norp_fit_coef,fit_freq);

c_list = matplotlib.colors.LinearSegmentedColormap.from_list(LAB.cm.datad['gist_rainbow'], colors=['r','g','b'], N = freqAmount);
saveFreqs = numpy.linspace(0,freqAmount-1,freqAmount,dtype=numpy.uint8)


# creating corrPlot

LAB.clf();
fig = LAB.figure(1,figsize = (20,16));
commonTitle = 'Correlation plot ' + startTime.strftime('%Y %B %d') + '. Stokes I, V at 4-8 GHz';
fig.suptitle(commonTitle,fontsize=8);

sub = fig.add_subplot(1,1,1);
sub.set_ylabel('correlation coefficient');
sub.set_xlabel('UT');
sub.xaxis.set_major_locator(LAB.MultipleLocator(dt_major));
sub.xaxis.set_major_formatter(LAB.FuncFormatter(hhmm_format));
sub.xaxis.set_minor_locator(LAB.MultipleLocator(dt_minor));
sub.set_xlim((0.*3600.,10.*3600.));
sub.set_ylim((-0.01,0.03));

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
LAB.savefig(corrPlotName_png);        
    

# writing srh_cp file
    
dataFormat = str(times.shape[1]) + 'D';
freqsColumn = fits.Column(name='frequencies',format='D',array=freqList[saveFreqs]);
timeColumn = fits.Column(name='time',format=dataFormat,array=saveTime);
IColumn = fits.Column(name='I',format=dataFormat,array=saveI);
VColumn = fits.Column(name='V',format=dataFormat,array=saveV);

fTableHdu = fits.BinTableHDU.from_columns([freqsColumn]);
dTableHdu = fits.BinTableHDU.from_columns([timeColumn, IColumn, VColumn]);
pHeader = fits.Header();
pHeader['DATE-OBS']     = sF.hduList[0].header['DATE-OBS'];
pHeader['TIME-OBS']     = sF.hduList[0].header['TIME-OBS'];
pHeader['INSTRUME']     = sF.hduList[0].header['INSTRUME'];
pHeader['ORIGIN']       = sF.hduList[0].header['ORIGIN'];
pHeader['OBS-LAT']      = sF.hduList[0].header['OBS-LAT'];
pHeader['OBS-LONG']     = sF.hduList[0].header['OBS-LONG'];
pHeader['OBS-ALT']      = sF.hduList[0].header['OBS-ALT'];
pHeader['FR_CHAN']      = sF.hduList[0].header['FR_CHAN'];
pHdu = fits.PrimaryHDU(header=pHeader);
hduList = fits.HDUList([pHdu, fTableHdu, dTableHdu]);
hduList.writeto('srh_cp_' + dateName + '.fits',clobber=True);
hduList.close();


# writing srh_amp file

dataFormat = str(times.shape[1]*48) + 'D';
rAmp_r = numpy.reshape(rAmp, (rAmp.shape[0], rAmp.shape[1]*rAmp.shape[2]))
lAmp_r = numpy.reshape(lAmp, (lAmp.shape[0], lAmp.shape[1]*lAmp.shape[2]))
rAmpColumn = fits.Column(name='AMP_RCP',format=dataFormat,array=rAmp_r);
lAmpColumn = fits.Column(name='AMP_LCP',format=dataFormat,array=lAmp_r);

antColumn = fits.Column(name='antenna',format='E',array=ant);

tTableHdu = fits.BinTableHDU.from_columns([timeColumn])
aTableHdu = fits.BinTableHDU.from_columns([antColumn])
dTableHdu = fits.BinTableHDU.from_columns([rAmpColumn,lAmpColumn]);


hduList = fits.HDUList([pHdu, fTableHdu, tTableHdu, aTableHdu, dTableHdu]);
hduList.writeto('srh_amp_' + dateName + '.fits',clobber=True);
hduList.close();


# writing srh_redVis file

dataFormat = str(times.shape[1]*rVis.shape[2]) + 'C';
rVis_r = numpy.reshape(rVis, (rVis.shape[0], rVis.shape[1]*rVis.shape[2]))
lVis_r = numpy.reshape(lVis, (lVis.shape[0], lVis.shape[1]*lVis.shape[2]))
rVisColumn = fits.Column(name='VIS_RCP',format=dataFormat,array=rVis_r);
lVisColumn = fits.Column(name='VIS_LCP',format=dataFormat,array=lVis_r);

antAColumn = fits.Column(name='antennaA',format='E',array=antA);
antBColumn = fits.Column(name='antennaB',format='E',array=antB);

tTableHdu = fits.BinTableHDU.from_columns([timeColumn])
aTableHdu = fits.BinTableHDU.from_columns([antAColumn, antBColumn])
dTableHdu = fits.BinTableHDU.from_columns([rVisColumn,lVisColumn]);

hduList = fits.HDUList([pHdu, fTableHdu, tTableHdu, aTableHdu, dTableHdu]);
hduList.writeto('srh_redVis_' + dateName + '.fits',clobber=True);
hduList.close();


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
#fdSrhImageName = open('srhImageName.txt','w');
#fdSrhImageName.write(srhImageName_png);
#fdSrhImageName.close();
#fi = open('srhImageName.txt','rb');
#fd.storbinary('STOR /Public/sergey/corrPlots/srhImageName.txt',fi);
#fi.close();
#fd.quit();

#
#srhSavNameIDL = 'srh_'+ dateName + '.sav'
#srhSav = readsav(srhSavNameIDL)
#savImage = srhSav['showlcprcp'][::-1]
#savImageFormat = str(savImage.shape[1]) + 'E'
#savFreqsColumn = fits.Column(name='frequencies',format='E',array=srhSav['sav_frequencies'])
#savTimeColumn = fits.Column(name='time',format='E',array=srhSav['sav_time'])
#savImageColumn = fits.Column(name='LCP_RCP',format=savImageFormat,array=savImage)
#
#savFreqTableHdu = fits.BinTableHDU.from_columns([savFreqsColumn])
#savTimeTableHdu = fits.BinTableHDU.from_columns([savTimeColumn])
#savImageTableHdu = fits.BinTableHDU.from_columns([savImageColumn]);
#
#pHeader = fits.Header();
#pHeader['DATE-OBS']     = sF.hduList[0].header['DATE-OBS'];
#pHeader['TIME-OBS']     = sF.hduList[0].header['TIME-OBS'];
#pHeader['INSTRUME']     = sF.hduList[0].header['INSTRUME'];
#pHeader['ORIGIN']       = sF.hduList[0].header['ORIGIN'];
#pHeader['OBS-LAT']      = sF.hduList[0].header['OBS-LAT'];
#pHeader['OBS-LONG']     = sF.hduList[0].header['OBS-LONG'];
#pHeader['OBS-ALT']      = sF.hduList[0].header['OBS-ALT'];
#pHeader['FR_CHAN']      = sF.hduList[0].header['FR_CHAN'];
#pHdu = fits.PrimaryHDU(header=pHeader);
#hduList = fits.HDUList([pHdu, savFreqTableHdu, savTimeTableHdu, savImageTableHdu]);
#hduList.writeto('srh_image_' + dateName + '.fits',clobber=True);
#hduList.close();

