#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 16:11:04 2020

@author: mariagloba
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 01:47:01 2016

@author: Sergey
"""
import os, fnmatch;
import srhFitsFile;
import numpy as NP;
import datetime as DT;
import pylab as PL;
import matplotlib;
from astropy.io import fits

dt_major = 3600.;
dt_minor = 900.;
freq_first = 2;
freq_intervals = 5;

def freq_format(f, pos):
  return '%3.1f' % ((f/FIT_FREQ_SIZE*4. + 4.));

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

def ihhmm_format(t, pos):
  t *= TIME_DELTA
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

os.environ['DISPLAY'] = ':0.0';

# Reading fits files

dateName = '20200421'
fitNames = findFits('/home/svlesovoi/SRH/2020/04/','mf_' + dateName + '*.fit')

fitNames.sort();
fitNames = fitNames[1:];
sF = srhFitsFile.SrhFitsFile(fitNames[0], 256)
for fName in fitNames[1:]:
    sF.append(fName)
freqAmount  = sF.freqListLength;
freqList = sF.freqList;
startTime= DT.datetime.strptime(sF.hduList[0].header['DATE-OBS'] + ' ' + sF.hduList[0].header['TIME-OBS'], '%Y/%m/%d %H:%M:%S.%f');
time = sF.freqTime;
fluxRcp = NP.mean(NP.abs(sF.visRcp[:,:,0:511]),axis=2);
fluxLcp = NP.mean(NP.abs(sF.visLcp[:,:,0:511]),axis=2);

vNorm = NP.mean(fluxRcp[:,:],axis=1) / NP.mean(fluxLcp[:,:],axis=1);

c_list = matplotlib.colors.LinearSegmentedColormap.from_list(PL.cm.datad['gist_rainbow'], colors=['r','g','b'], N = freqAmount);
saveFreqs = NP.linspace(0,freqAmount-1,freqAmount,dtype=NP.uint8)

PL.clf();
fig = PL.figure(1,figsize = (20,16));
commonTitle = 'Correlation plot ' + startTime.strftime('%Y %B %d') + '. Stokes I, V at 4-8 GHz';
fig.suptitle(commonTitle,fontsize=8);

sub = fig.add_subplot(1,1,1);
sub.set_ylabel('correlation coefficient');
sub.set_xlabel('UT');
sub.xaxis.set_major_locator(PL.MultipleLocator(dt_major));
sub.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format));
sub.xaxis.set_minor_locator(PL.MultipleLocator(dt_minor));
sub.set_xlim((0.*3600.,10.*3600.));
sub.set_ylim((-0.01,0.03));

times = NP.zeros((saveFreqs.shape[0],time.shape[1]));
i_map = NP.zeros((saveFreqs.shape[0],time.shape[1]));
v_map = NP.zeros((saveFreqs.shape[0],time.shape[1]));
for ff in saveFreqs:
    sub.scatter(time[ff,:],0.5*(fluxRcp[ff,:] + vNorm[ff]*fluxLcp[ff,:]),color=c_list(ff), s=0.5, linewidths = 0)
    sub.scatter(time[ff,:],0.5*(fluxRcp[ff,:] - vNorm[ff]*fluxLcp[ff,:]),color=c_list(ff), s=0.5, linewidths = 0)
    times[ff,:] = time[ff,:];
    i_map[ff,:] = 0.5*(fluxRcp[ff,:] + vNorm[ff]*fluxLcp[ff,:]);
    v_map[ff,:] = 0.5*(fluxRcp[ff,:] - vNorm[ff]*fluxLcp[ff,:]);
        
dataFormat = str(times.shape[1]) + 'D';
freqsColumn = fits.Column(name='frequencies',format='D',array=freqList[saveFreqs]);
timeColumn = fits.Column(name='time',format=dataFormat,array=times);
IColumn = fits.Column(name='I',format=dataFormat,array=i_map);
VColumn = fits.Column(name='V',format=dataFormat,array=v_map);

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
#hduList.writeto('srh_cp_' + dateName + '.fits',clobber=True);
hduList.close();

corrPlotName = 'fCorrPlot'+ dateName;
corrPlotName_png = corrPlotName + '.png';
PL.savefig(corrPlotName_png, dpi = 200);

# Making corrMap

TIME_DELTA = 10#8.97
FIT_FREQ_SIZE = 1500
FIT_TIME_SIZE = 3600

timeOfStart= DT.datetime.strptime(sF.hduList[0].header['DATE-OBS'], '%Y/%m/%d');
hdr = {'cdelt1':1, 'cdelt2':1}

freqs = freqList[saveFreqs]

NP.clip(i_map, -0.2, 0.2)
NP.clip(v_map, -0.2, 0.2)
src_t0 = times[0,0]
src_t1 = max(times[0,:])
cal_mean_zero = NP.zeros(i_map.shape[0])

fit_dt = 10
fit_t0 = 0 * fit_dt * 3600
fit_t1 = 10 * fit_dt * 3600
fit_time_size = i_map.shape[1]
fit_freqs = (4. + 4.*NP.arange(FIT_FREQ_SIZE)/FIT_FREQ_SIZE)*1e3
fit_times = (src_t0 +fit_dt* NP.arange(fit_time_size))

cal_mean = NP.mean(i_map)
cal_min = NP.min(i_map)

if (i_map.shape[1] > 1300):
    cal_min_column = i_map[:,1120].copy()
    cal_column = i_map[:,1200].copy()
else:
    cal_min_column = cal_mean_zero.copy()
    cal_column = NP.mean(i_map,1).copy()

fit_i_map = NP.zeros((FIT_FREQ_SIZE,FIT_TIME_SIZE))
fit_v_map = NP.zeros((FIT_FREQ_SIZE,FIT_TIME_SIZE))
for tt in range(fit_time_size):
    i_map[:,tt] -=cal_min_column
    i_map[:,tt] =i_map[:,tt] / cal_column * cal_mean * 3
    fit_tt = int(times[0, tt] / 10 + .5)
    fit_i_map[:,fit_tt] = NP.interp(fit_freqs,freqs,i_map[:,tt])
    fit_v_map[:,fit_tt] = NP.interp(fit_freqs,freqs,v_map[:,tt])

cal_mean_column= 0
result_time_delta = (times[0,10] - times[0,9])
dt_major = 360#3600/result_time_delta;
dt_minor = 72#900/result_time_delta;
df_major = 500/(4000/FIT_FREQ_SIZE)
df_minor = 50/(4000/FIT_FREQ_SIZE)

PL.clf();
fig = PL.figure(figsize=(16,12));
commonTitle = 'Summed cross-correlation amplitude for '  + timeOfStart.strftime('%Y %m %d')
fig.suptitle(commonTitle,fontsize='xx-large')

spI = fig.add_subplot(2,1,1,title='RCP + LCP')
spI.set_ylabel('Frequency [GHz]')
spI.xaxis.set_major_locator(PL.MultipleLocator(dt_major));
spI.xaxis.set_minor_locator(PL.MultipleLocator(dt_minor));
spI.xaxis.set_major_formatter(PL.FuncFormatter(ihhmm_format))
spI.yaxis.set_major_locator(PL.MultipleLocator(df_major));
spI.yaxis.set_minor_locator(PL.MultipleLocator(df_minor));
spI.yaxis.set_major_formatter(PL.FuncFormatter(freq_format))
spI.imshow(fit_i_map*300, origin='lower',cmap='ocean',vmax=2.0)
fig.colorbar(spI.images[0],label='Amplitude [arb. units]');

spV = fig.add_subplot(2,1,2,title='RCP - LCP')
spV.set_ylabel('Frequency [GHz]')
spV.xaxis.set_major_locator(PL.MultipleLocator(dt_major));
spV.xaxis.set_minor_locator(PL.MultipleLocator(dt_minor));
spV.xaxis.set_major_formatter(PL.FuncFormatter(ihhmm_format))
spV.yaxis.set_major_formatter(PL.FuncFormatter(freq_format))
spV.imshow(fit_v_map*300, origin='lower',cmap='ocean',vmin=-.5,vmax=.5)
spV.set_xlabel('Time [UT]', horizontalalignment='right')
spV.xaxis.set_label_coords(1,-0.1)
fig.colorbar(spV.images[0],label='Amplitude [arb. units]');
