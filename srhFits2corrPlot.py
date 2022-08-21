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


os.environ['DISPLAY'] = ':0.0';

dateName = DT.datetime.now().strftime("%Y%m%d");
fitNames = findFits('/home/svlesovoi/SRH/2016/07/23/','mf_*.fit')

fitNames.sort();
fitNames = fitNames[1:];
#fdSrhImageName = open('/home/serg/SRH/mf_currentFileName.txt','w');
#fdSrhImageName.write(fitNames[-1]);
#fdSrhImageName.close();
visFirst = 0
visLast = 512
saveFreqs = numpy.array([2,5,8,11,14])
scaleII = 2.55e-3
ampScale = False

for fName in fitNames:
    sF = srhFitsFile.SrhFitsFile(fName, 256)
    if (fName == fitNames[0]):
        freqAmount  = sF.freqListLength;
        freqList = sF.freqList;
        startTime= DT.datetime.strptime(sF.hduList[0].header['DATE-OBS'] + ' ' + sF.hduList[0].header['TIME-OBS'], '%Y/%m/%d %H:%M:%S.%f');
        times = sF.freqTime;
        r_corr = numpy.sin(numpy.pi/2.*numpy.abs(sF.visRcp[saveFreqs,:,visFirst:visLast]))
        r_amp = sF.ampRcp[saveFreqs,:,:]
        if (ampScale):
            for i in range(saveFreqs.shape[0]):
                for j in range(16):
                    for k in range(32):
                        II = numpy.sqrt(r_amp[i,:,j]) * numpy.sqrt(r_amp[i,:,k + 16]) * scaleII
                        r_corr[i,:,j*k] *= II
        corrRcp = numpy.mean(r_corr,axis=2);
#        ampRcp = numpy.mean(r_amp,axis=2)
        ampRcp = r_amp
        r_corr = numpy.sin(numpy.pi/2.*numpy.abs(sF.visLcp[saveFreqs,:,visFirst:visLast]))
        r_amp = sF.ampLcp[saveFreqs,:,:]
        corr_0_192_49 = r_corr[0,:,0]
        if (ampScale):
            for i in range(saveFreqs.shape[0]):
                for j in range(16):
                    for k in range(32):
                        II = numpy.sqrt(r_amp[i,:,j]) * numpy.sqrt(r_amp[i,:,k + 16]) * scaleII
                        r_corr[i,:,j*k] *= II
        corrLcp = numpy.mean(r_corr,axis=2);
        #ampLcp = numpy.mean(r_amp,axis=2)
        ampLcp = r_amp
        r_corr = 0;
        r_amp = 0;
    else:
        if sF.freqTime[0,1] - times[0,-1] > 600.:
            insTime = numpy.zeros((freqAmount,2),dtype=float);
            insRcp = numpy.zeros((freqAmount,2),dtype=float);
            insLcp = numpy.zeros((freqAmount,2),dtype=float);
            for fff in range(freqAmount):
                insTime[fff,0] = times[fff,-1];
                insTime[fff,1] = sF.freqTime[fff,0]
            times = numpy.concatenate((times,insTime),axis=1);
            corrRcp = numpy.concatenate((corrRcp,insRcp),axis=1);
            corrLcp = numpy.concatenate((corrLcp,insLcp),axis=1);
        times = numpy.concatenate((times,sF.freqTime),axis=1);
        r_corr = numpy.sin(numpy.pi/2.*numpy.abs(sF.visRcp[saveFreqs,:,visFirst:visLast]))
        r_amp = sF.ampRcp[saveFreqs,:,:]
        if (ampScale):
            for i in range(saveFreqs.shape[0]):
                for j in range(16):
                    for k in range(32):
                        II = numpy.sqrt(r_amp[i,:,j]) * numpy.sqrt(r_amp[i,:,k + 16]) * scaleII
                        r_corr[i,:,j*k] *= II
        corrRcp = numpy.concatenate((corrRcp, numpy.mean(r_corr,axis=2)),axis=1);
#        ampRcp = numpy.concatenate((ampRcp, numpy.mean(r_amp,axis=2)),axis=1)
        ampRcp = numpy.concatenate((ampRcp, r_amp),axis=1)
        corr_0_192_49 = numpy.concatenate((corr_0_192_49, r_corr[0,:,0]))

        r_corr = numpy.sin(numpy.pi/2.*numpy.abs(sF.visLcp[saveFreqs,:,visFirst:visLast]))
        r_amp = sF.ampLcp[saveFreqs,:,:]
        if (ampScale):
            for i in range(saveFreqs.shape[0]):
                for j in range(16):
                    for k in range(32):
                        II = numpy.sqrt(r_amp[i,:,j]) * numpy.sqrt(r_amp[i,:,k + 16]) * scaleII
                        r_corr[i,:,j*k] *= II
        corrLcp = numpy.concatenate((corrLcp, numpy.mean(r_corr,axis=2)),axis=1);
#        ampLcp = numpy.concatenate((ampLcp, numpy.mean(r_amp,axis=2)),axis=1)
        ampLcp = numpy.concatenate((ampLcp, r_amp),axis=1)
        r_corr = 0;
        r_amp = 0;
    sF.close();

vNorm = numpy.mean(corrRcp[:,:],axis=1) / numpy.mean(corrLcp[:,:],axis=1);

c_list = matplotlib.colors.LinearSegmentedColormap.from_list(LAB.cm.datad['gist_rainbow'], colors=['r','g','b'], N = freqAmount);
saveFreqs = numpy.array([2,5,8,11,14])

#LAB.clf();
fig = LAB.figure(figsize = (20,16));
commonTitle = 'Correlation plot ' + startTime.strftime('%Y %B %d') + '. Stokes I, V at 4-8 GHz';
#for ff in saveFreqs:
#    commonTitle += "%1.1f" % (freqList[ff]*1e-3);
#    if (ff < freqAmount - 1):
#        commonTitle += ", ";
#commonTitle += " GHz";
fig.suptitle(commonTitle,fontsize=8);

sub = fig.add_subplot(1,1,1);
sub.set_ylabel('correlation coefficient');
sub.set_xlabel('UT');
sub.xaxis.set_major_locator(LAB.MultipleLocator(dt_major));
sub.xaxis.set_major_formatter(LAB.FuncFormatter(hhmm_format));
sub.xaxis.set_minor_locator(LAB.MultipleLocator(dt_minor));
sub.set_xlim((0.*3600.,10.*3600.));
sub.set_ylim((-0.1,0.60));

saveTime = numpy.zeros((saveFreqs.shape[0],times.shape[1]));
saveI = numpy.zeros((saveFreqs.shape[0],times.shape[1]));
saveV = numpy.zeros((saveFreqs.shape[0],times.shape[1]));
#for ff in saveFreqs:
for i in range(saveFreqs.shape[0]):
    ff = saveFreqs[i]
    sub.plot(times[ff,:],0.5*(corrRcp[i,:] + vNorm[i]*corrLcp[i,:]),color=c_list(ff))#,ls='None',marker='.',markersize=.3)
    sub.plot(times[ff,:],0.5*(corrRcp[i,:] - vNorm[i]*corrLcp[i,:]),color=c_list(ff))#,ls='None',marker='.',markersize=.3)
    saveTime[i,:] = times[ff,:];
    saveI[i,:] = 0.5*(corrRcp[i,:] + vNorm[i]*corrLcp[i,:]);
    saveV[i,:] = 0.5*(corrRcp[i,:] - vNorm[i]*corrLcp[i,:]);
        

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

corrPlotName = 'fCorrPlot'+ dateName;
corrPlotName_png = corrPlotName + '.png';
LAB.savefig(corrPlotName_png);

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
#
#fd.quit();
