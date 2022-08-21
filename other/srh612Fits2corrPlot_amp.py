# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 01:47:01 2016

@author: Sergey
"""
import os, fnmatch;
from srhFitsFile612 import SrhFitsFile
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

def freq2flux(freq):
    return (freq - 2800000)*40/3000000 + 80
    
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
dateName = '20211210'
#fitNames = findFits('/home/sergeyvlesovoi/NFS_TEST/SRH/SRH0612/' + dateName ,'*.fit')  
#fitNames = findFits('/home/sergeyvlesovoi/SRH_0612/' + dateName ,'*.fit')  
#fitNames = findFits('SRH0612/' + dateName ,'*.fit')  
fitNames = findFits('/home/sergey_lesovoi/SRH_DATA/SRH/SRH0612/' + dateName ,'*.fit')  
fitNames.sort();
fitNames = fitNames[1:];

visFirst = 32*0
visLast = 32*16 - 1

for fName in fitNames:
    print(fName)
    sF = SrhFitsFile(fName, 256)
    if (fName == fitNames[0]):
        freqAmount  = sF.freqListLength;
        freqList = sF.freqList;
        startTime= DT.datetime.strptime(sF.hduList[0].header['DATE-OBS'] + ' ' + sF.hduList[0].header['TIME-OBS'], '%Y-%m-%d %H:%M:%S');
        times = sF.freqTime;
        r_real = abs(sF.visRcp.real[:,:,0:8192])
        r_imag = abs(sF.visRcp.imag[:,:,0:8192])
        corrFluxRcp = numpy.mean(numpy.sqrt(r_real**2. + r_imag**2.),axis=2);
        r_real = abs(sF.visLcp.real[:,:,0:3007])
        r_imag = abs(sF.visLcp.imag[:,:,0:3007])
        corrFluxLcp = numpy.mean(numpy.sqrt(r_real**2. + r_imag**2.),axis=2);
        r_real = 0;
        r_imag = 0;
        ampFluxRcp = numpy.mean(sF.ampRcp, axis = 2)
        ampFluxLcp = numpy.mean(sF.ampLcp, axis = 2)
        
        
    else:
        if sF.freqTime[0,1] - times[0,-1] > 600.:
            insTime = numpy.zeros((freqAmount,2),dtype=float);
            insRcp = numpy.zeros((freqAmount,2),dtype=float);
            insLcp = numpy.zeros((freqAmount,2),dtype=float);
            for fff in range(freqAmount):
                insTime[fff,0] = times[fff,-1];
                insTime[fff,1] = sF.freqTime[fff,0]
            times = numpy.concatenate((times,insTime),axis=1);
            corrFluxRcp = numpy.concatenate((corrFluxRcp,insRcp),axis=1);
            corrFluxLcp = numpy.concatenate((corrFluxLcp,insLcp),axis=1);
        times = numpy.concatenate((times,sF.freqTime),axis=1);
        r_real = abs(sF.visRcp.real[:,:,0:3007])
        r_imag = abs(sF.visRcp.imag[:,:,0:3007])
        corrFluxRcp = numpy.concatenate((corrFluxRcp, numpy.mean(numpy.sqrt(r_real**2. + r_imag**2.),axis=2)),axis=1);
        r_real = abs(sF.visLcp.real[:,:,0:3007])
        r_imag = abs(sF.visLcp.imag[:,:,0:3007])
        corrFluxLcp = numpy.concatenate((corrFluxLcp, numpy.mean(numpy.sqrt(r_real**2. + r_imag**2.),axis=2)),axis=1);
        r_real = 0;
        r_imag = 0;
        ampFluxRcp = numpy.concatenate((ampFluxRcp, numpy.mean(sF.ampRcp, axis = 2)), axis = 1)
        ampFluxLcp = numpy.concatenate((ampFluxLcp, numpy.mean(sF.ampLcp, axis = 2)), axis = 1)
        
    sF.close();

vNorm = numpy.mean(corrFluxRcp[:,:],axis=1) / numpy.mean(corrFluxLcp[:,:],axis=1);
vNorm_amp = numpy.mean(ampFluxRcp[:,:],axis=1) / numpy.mean(ampFluxLcp[:,:],axis=1);

c_list = matplotlib.colors.LinearSegmentedColormap.from_list(LAB.cm.datad['gist_rainbow'], colors=['r','g','b'], N = freqAmount);
saveFreqs = numpy.linspace(0,freqAmount-1,freqAmount,dtype=numpy.uint8)

# creating corrPlot

fig = LAB.figure(figsize = (16,8));
fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.05)
fig.tight_layout()
commonTitle = 'Correlation plot ' + startTime.strftime('%Y %B %d') + '. Stokes I, V at 6-12 GHz';
fig.suptitle(commonTitle,fontsize=14);

sub = fig.add_subplot(1,1,1);
sub.set_ylabel('correlation coefficient');
sub.set_xlabel('UT');
sub.xaxis.set_major_locator(LAB.MultipleLocator(dt_major));
sub.xaxis.set_major_formatter(LAB.FuncFormatter(hhmm_format));
sub.xaxis.set_minor_locator(LAB.MultipleLocator(dt_minor));
sub.set_xlim((0.*3600.,11.*3600.));
sub.set_ylim((-0.002,0.02));

saveTime = numpy.zeros((saveFreqs.shape[0],times.shape[1]));
saveI = numpy.zeros((saveFreqs.shape[0],times.shape[1]));
saveV = numpy.zeros((saveFreqs.shape[0],times.shape[1]));
zeroLevel = 0.0009
for ff in saveFreqs:
    sub.scatter(times[ff,:],(0.5*(corrFluxRcp[ff,:] + vNorm[ff]*corrFluxLcp[ff,:])-zeroLevel),color=c_list(ff), s=0.5, linewidths = 0,label='%d MHz'%(freqList[ff]*1e-3))
    sub.scatter(times[ff,:],0.5*(corrFluxRcp[ff,:] - vNorm[ff]*corrFluxLcp[ff,:]),color=c_list(ff), s=0.5, linewidths = 0)
    sub.legend(markerscale=10)
    saveTime[ff,:] = times[ff,:];
    saveI[ff,:] = 0.5*(corrFluxRcp[ff,:] + vNorm[ff]*corrFluxLcp[ff,:]);
    saveV[ff,:] = 0.5*(corrFluxRcp[ff,:] - vNorm[ff]*corrFluxLcp[ff,:]);
corrPlotName = 'fCorrPlot'+ dateName;
corrPlotName_png = corrPlotName + '.png';
LAB.savefig(corrPlotName_png); 

# creating fluxPlot

fig = LAB.figure(figsize = (16,8));
fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
fig.tight_layout()
commonTitle = 'Flux plot ' + startTime.strftime('%Y %B %d') + '. Stokes I, V at 3-6 GHz';
fig.suptitle(commonTitle,fontsize=14);

sub = fig.add_subplot(1,1,1);
sub.set_ylabel('flux, a.u.');
sub.set_xlabel('UT');
sub.xaxis.set_major_locator(LAB.MultipleLocator(dt_major));
sub.xaxis.set_major_formatter(LAB.FuncFormatter(hhmm_format));
sub.xaxis.set_minor_locator(LAB.MultipleLocator(dt_minor));
sub.set_xlim((0.*3600.,11.*3600.));
#sub.set_ylim((-2e3, 1e4));
sub.set_ylim((-2e1, 3e2));

saveI_amp = numpy.zeros((saveFreqs.shape[0],ampFluxRcp.shape[1]));
saveV_amp = numpy.zeros((saveFreqs.shape[0],ampFluxRcp.shape[1]));
ampLength = ampFluxRcp[ff,:].shape[0]
for ff in saveFreqs:
    sub.scatter(times[ff,0:ampLength],(0.5*(ampFluxRcp[ff,:] + vNorm_amp[ff]*ampFluxLcp[ff,:]))*freq2flux(freqList[ff]),color=c_list(ff), s=1.5, linewidths = 0,label='%d MHz'%(freqList[ff]*1e-3))
    sub.scatter(times[ff,0:ampLength],(0.5*(ampFluxRcp[ff,:] - vNorm_amp[ff]*ampFluxLcp[ff,:]))*freq2flux(freqList[ff]),color=c_list(ff), s=1.5, linewidths = 0)
    sub.legend(markerscale=5)
    saveI_amp[ff,:] = 0.5*(ampFluxRcp[ff,:] + vNorm_amp[ff]*ampFluxLcp[ff,:]);
    saveV_amp[ff,:] = 0.5*(ampFluxRcp[ff,:] - vNorm_amp[ff]*ampFluxLcp[ff,:]);
fluxPlotName = 'fFluxPlot'+ dateName;
fluxPlotName_png = fluxPlotName + '.png';
LAB.savefig(fluxPlotName_png);        
    

# writing srh_cp file
    
dataFormat = str(times.shape[1]) + 'D';
freqsColumn = fits.Column(name='frequencies',format='D',array=freqList[saveFreqs]);
timeColumn = fits.Column(name='time',format=dataFormat,array=saveTime);
IColumn = fits.Column(name='I',format=dataFormat,array=saveI);
VColumn = fits.Column(name='V',format=dataFormat,array=saveV);
IAmpColumn = fits.Column(name='flux_I',format=dataFormat,array=saveI_amp);
VAmpColumn = fits.Column(name='flux_V',format=dataFormat,array=saveV_amp);

fTableHdu = fits.BinTableHDU.from_columns([freqsColumn]);
dTableHdu = fits.BinTableHDU.from_columns([timeColumn, IColumn, VColumn, IAmpColumn, VAmpColumn]);
pHeader = fits.Header();
pHeader['DATE-OBS']     = sF.hduList[0].header['DATE-OBS'];
pHeader['TIME-OBS']     = sF.hduList[0].header['TIME-OBS'];
pHeader['INSTRUME']     = 'SRH'
pHeader['ORIGIN']       = 'ISTP'
pHeader['OBS-LAT']      = '51.759'
pHeader['OBS-LONG']     = '102.217'
pHeader['OBS-ALT']      = '799'
pHeader['FR_CHAN']      = '10'
pHdu = fits.PrimaryHDU(header=pHeader);
hduList = fits.HDUList([pHdu, fTableHdu, dTableHdu]);
hduList.writeto('srh_cp_' + dateName + '.fits',clobber=True);
hduList.close();

