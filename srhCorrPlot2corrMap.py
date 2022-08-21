#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 06:58:08 2017

@author: sergey
"""

from astropy.io import fits
import pylab as PL
import numpy as NP
import datetime as DT
from scipy import stats
import ftplib;

#TIME_DELTA = 16.13
TIME_DELTA = 10#8.97
FIT_FREQ_SIZE = 1500
FIT_TIME_SIZE = 3600

def freq_format(f, pos):
  return '%3.1f' % ((f/FIT_FREQ_SIZE*4. + 4.));

def ihhmm_format(t, pos):
  t *= TIME_DELTA
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

dateName = DT.datetime.now().strftime("%Y%m%d");
cpHdu = fits.open('/home/serg/py/srh_cp_'+ dateName + '.fits')
#dateName = '20191120';
#cpHdu = fits.open('/home/svlesovoi/SRH/lightcurve/2019/srh_cp_'+ dateName + '.fits')
timeOfStart= DT.datetime.strptime(cpHdu[0].header['DATE-OBS'], '%Y/%m/%d');
hdr = {'cdelt1':1, 'cdelt2':1}

freqs = cpHdu[1].data['frequencies']
times = cpHdu[2].data['time']
i_map = cpHdu[2].data['I']
v_map = cpHdu[2].data['V']
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
    i_map[:,tt] =i_map[:,tt] / cal_column * cal_mean
    fit_tt = int(times[0, tt] / 10 + .5)
    fit_i_map[:,fit_tt] = NP.interp(fit_freqs,freqs,i_map[:,tt])
    fit_v_map[:,fit_tt] = NP.interp(fit_freqs,freqs,v_map[:,tt])

cal_mean_column= 0
result_time_delta = (times[0,10] - times[0,9])
dt_major = 360#3600/result_time_delta;
dt_minor = 72#900/result_time_delta;
df_major = 500/(4000/FIT_FREQ_SIZE)
df_minor = 50/(4000/FIT_FREQ_SIZE)

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

fig.savefig('srh_corrmap_'+ dateName + '.png');

fd = ftplib.FTP('192.168.0.1','sergey','jill21');
fi = open('srh_corrmap_' + dateName + '.png','rb');
fd.storbinary('STOR /Public/sergey/corrPlots/' + 'srh_corrmap_' + dateName + '.png',fi);
fi.close();

corrMapName = 'srh_corrmap_'+ dateName;
corrMapName_png = corrMapName + '.png';

fdCPlotName = open('corrMapName.txt','w');
fdCPlotName.write(corrMapName_png);
fdCPlotName.close();

fi = open('corrMapName.txt','rb');
fd.storbinary('STOR /Public/sergey/corrPlots/corrMapName.txt',fi);
fi.close();
