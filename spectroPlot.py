#!/usr/bin/python
# -*- coding: utf-8 -*-
import pylab as LAB;
import pyfits as FITS;
import datetime as DT;
import matplotlib;
from matplotlib.ticker import MultipleLocator, FormatStrFormatter;
import os, fnmatch;
import numpy;
from scipy import *;
import ftplib;
import string;
from optparse import OptionParser;
import matplotlib.font_manager as fm;
import sunpy.image as SUNIM;

def freq_format(f, pos):
  return '%3.1f' % ((f/400.*22. + 2.));

def hhmmss_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  t -= mm*60.;
  ss = (int)(t);
  return '%02d:%02d:%02d' % (hh,mm,ss);

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

def ihhmm_format(t, pos):
  t *= 60.;
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

def findFits(path, pattern, size):
  result = [];
  for root, dirs, files in os.walk(path):
    for basename in files:
      if fnmatch.fnmatch(basename,pattern):
	if os.path.getsize(os.path.join(root,basename)) > size:
	  result.append(os.path.join(root,basename));
  return result;

os.environ['DISPLAY'] = ':0.0';

parser = OptionParser();
parser.add_option("-d", "--date", dest="inDate",   default="20110225");
parser.add_option("-s", "--size", dest="fitsSize", type="int",  default=2880);
#parser.add_option("-s", "--size", dest="fitsSize", type="int",  default=37440);
(cl_options,cl_args) = parser.parse_args();
dateName = cl_options.inDate;

fdName = open('/home/sergey/SSRT/data/spectro/aCurrentData.txt','rb');
#fdName = open('/home/serg/SSRT/data/spectro/aCurrentData.txt','rb');
dateName = fdName.readline();
fdName.close();

fitName = findFits('/home/sergey/SSRT/data/spectro/','*smf*' + dateName + '*.fit',cl_options.fitsSize);
#fitName = findFits('/home/serg/SSRT/data/spectro/','*smf*' + dateName + '*.fit',cl_options.fitsSize);

fitName.sort();

i = 0;
for fName in fitName:
  hdu = FITS.open(fName);

  if i == 0:
    timeOfStart= DT.datetime.strptime(hdu[0].header['DATE-OBS'] + ' ' + hdu[0].header['TIME-OBS'], '%Y/%m/%d %H:%M:%S.%f');
    freqs	= hdu[1].data.field('FREQUENCY')*1e-9;
    freq_time	= hdu[1].data.field('TIME');
    rcp_data	= hdu[1].data.field('RCP');
    lcp_data	= hdu[1].data.field('LCP');
    freqAmount  = freqs.shape[0];
  else:
    freq_time = concatenate(( freq_time ,hdu[1].data.field('TIME')),axis=1);
    rcp_data  = concatenate((rcp_data,hdu[1].data.field('RCP')),axis=1);
    lcp_data  = concatenate((lcp_data,hdu[1].data.field('LCP')),axis=1);
  i = i + 1;
  hdu.close();

i = 0;

L = rcp_data.shape[1];
F = rcp_data.shape[0];

fdNorpData = open('norp_data.txt','rb');
norpLines = fdNorpData.readlines();
norpData = [];
for i in norpLines:
	norpData.append([float(j) for j in i.split()]);
norpData = array(norpData);
fdNorpData.close();

norp_freq = norpData[:,0];
norp_flux = norpData[:,1];
norp_fit_coef = polyfit(norp_freq,norp_flux,2);

fit_freq = 25.*arange(100.)/100.;
fit_flux = polyval(norp_fit_coef,fit_freq);

sp_2_24_flux = zeros(F);
for i in range(0,F):
	sp_2_24_flux[i] = fit_flux[where(abs(fit_freq - freqs[i]) < 0.125)];

clbSky = [];
fdClbSky = open('/home/sergey/py/sclb_sky.txt','rb');
#fdClbSky = open('/home/serg/py/sclb_sky.txt','rb');
for clbSkyLine in fdClbSky.readlines():
  clbSkyRow = clbSkyLine.split();
  for clbSkyWord in clbSkyRow:
    clbSky.append(float(clbSkyWord));
min_sky = array(clbSky).reshape(F,2).transpose();
fdClbSky.close();

clbSun = [];
fdClbSun = open('/home/sergey/py/sclb_sun.txt','rb');
#fdClbSun = open('/home/serg/py/sclb_sun.txt','rb');
for clbSunLine in fdClbSun.readlines():
  clbSunRow = clbSunLine.split();
  for clbSunWord in clbSunRow:
    clbSun.append(float(clbSunWord));
mean_sun = array(clbSun).reshape(F,2).transpose();
fdClbSun.close();

for i in range(0,F):
  lcp_data[i,:] = lcp_data[i,:] - min_sky[0,i];
  lcp_data[i,:] = lcp_data[i,:]/mean_sun[0,i] * sp_2_24_flux[i];
  rcp_data[i,:] = rcp_data[i,:] - min_sky[1,i];
  rcp_data[i,:] = rcp_data[i,:]/mean_sun[1,i] * sp_2_24_flux[i];

for i in range(1,L):
  if (freq_time[0,i] == -1.):
    freq_time[0,i] = freq_time[0,i - 1];

i_data = (rcp_data + lcp_data)/2.;
v_data = (rcp_data - lcp_data)/2.;

c_list = matplotlib.colors.LinearSegmentedColormap.from_list(LAB.cm.datad['gist_rainbow'], colors=['r','g','b'], N = F);

fig1 = LAB.figure(1,figsize=(20,12));
commonTitle = 'Siberian spectropolarimeter 2-24 GHz, ' + timeOfStart.strftime('%Y %m %d');
fig1.suptitle(commonTitle,fontsize='large');

dt_major = 3600;
dt_minor = 900;
t0 = 0.* 3600.;
t1 = 10.5* 3600.;
it0 = freq_time[0,0];
it1 = max(freq_time[0]);

zero_ind1 = where(abs(freq_time[0,:] - roll(freq_time[0,:],1)) > 10);
zero_ind1 = array(zero_ind1[0]);
zero_ind2 = where(abs(freq_time[0,:] - roll(freq_time[0,:],-1)) > 10);
zero_ind2 = array(zero_ind2[0]);
if (zero_ind1.shape > 0):
  i_data[zero_ind1,:] = 0.;
  v_data[zero_ind1,:] = 0.;
  
if (zero_ind2.shape > 0):
  i_data[zero_ind2,:] = 0.;
  v_data[zero_ind2,:] = 0.;

sp1 = fig1.add_subplot(2,1,1);
sp1.set_ylabel('RCP + LCP [s.f.u.]');
sp1.set_xlabel('UT');
sp1.xaxis.set_major_locator(MultipleLocator(dt_major));
sp1.xaxis.set_major_formatter(LAB.FuncFormatter(hhmm_format));
sp1.xaxis.set_minor_locator(MultipleLocator(dt_minor));
for i in range(0,F):
  sp1.plot(freq_time[0,:],i_data[i,:],color = c_list(float(i)/F));
sp1.set_xlim((t0,t1));
sp1.set_ylim((0.,3000.));
legendList = ["%1.2f" % freqs[0]];
for i in range(1,F):
  legendList.append("%1.2f" % freqs[i]);
prop = fm.FontProperties(size=10);
sp1.legend(legendList,prop=prop,ncol=2,title="GHz");

sp2 = fig1.add_subplot(2,1,2);
sp2.set_ylabel('RCP - LCP [s.f.u.]');
sp2.set_xlabel('UT');
sp2.xaxis.set_major_locator(MultipleLocator(dt_major));
sp2.xaxis.set_major_formatter(LAB.FuncFormatter(hhmm_format));
sp2.xaxis.set_minor_locator(MultipleLocator(dt_minor));
for i in range(0,F):
  sp2.plot(freq_time[0,:],v_data[i,:],color = c_list(float(i)/F));
sp2.set_xlim((t0,t1));
sp2.set_ylim((-1500.,1500.));

new_width  = int((t1 - t0)/60.);
new_height = 400;
new_t0 = int(it0 / 60.);
new_t1 = int(it1 / 60.);
fit_freq = (2. + 22.*arange(new_height)/new_height);

res_i = SUNIM.resample(i_data,(i_data.shape[0],new_t1 - new_t0));
res_v = SUNIM.resample(v_data,(v_data.shape[0],new_t1 - new_t0));
ii = zeros((new_height,new_width));
iv = zeros((new_height,new_width));

for i in range(new_t1 - new_t0):
	_flux = res_i[:,i];
	ii[:,i + new_t0] = interp(fit_freq,freqs,_flux);
	_flux = res_v[:,i];
	iv[:,i + new_t0] = interp(fit_freq,freqs,_flux);

fig2 = LAB.figure(2,figsize=(20,12));
commonTitle = 'Siberian spectropolarimeter 2-24 GHz, ' + timeOfStart.strftime('%Y %m %d');
fig2.suptitle(commonTitle,fontsize='large');

sp3 = fig2.add_subplot(2,1,1,title='RCP + LCP');
sp3.set_xlabel('UT');
sp3.set_ylabel('GHz');
sp3.xaxis.set_major_locator(MultipleLocator(dt_major/60.));
sp3.xaxis.set_major_formatter(LAB.FuncFormatter(ihhmm_format));
sp3.yaxis.set_major_formatter(LAB.FuncFormatter(freq_format));
sp3.yaxis.set_major_locator(MultipleLocator(new_height/11.));
sp3.imshow(ii,origin='lower',aspect=0.5);
fig2.colorbar(sp3.images[0]);

sp4 = fig2.add_subplot(2,1,2,title='RCP - LCP');
sp4.set_xlabel('UT');
sp4.set_ylabel('GHz');
sp4.xaxis.set_major_locator(MultipleLocator(dt_major/60.));
sp4.xaxis.set_major_formatter(LAB.FuncFormatter(ihhmm_format));
sp4.yaxis.set_major_formatter(LAB.FuncFormatter(freq_format));
sp4.yaxis.set_major_locator(MultipleLocator(new_height/11.));
sp4.imshow(iv,origin='lower',aspect=0.5);
fig2.colorbar(sp4.images[0]);

spectroPlotName = 'spectroPlot'+ DT.datetime.now().strftime("%Y%m%d");
spectroPlotName1_png = spectroPlotName + '_1.png';
spectroPlotName2_png = spectroPlotName + '_2.png';
spectroPlotName_txt = spectroPlotName + '.txt';

fig1.savefig(spectroPlotName1_png);
fig2.savefig(spectroPlotName2_png);

#LAB.show();

fdCPlotName = open('spectroPlotName.txt','wb');
fdCPlotName.write(spectroPlotName1_png);
fdCPlotName.write('\n');
fdCPlotName.write(spectroPlotName2_png);
fdCPlotName.write('\n');
fdCPlotName.close();

fd = ftplib.FTP('192.168.0.1','sergey','jill21');

fi = open(spectroPlotName1_png,'rb');
fd.storbinary('STOR /Public/sergey/spectroPlots/' + spectroPlotName1_png,fi);
fi.close();

fi = open(spectroPlotName2_png,'rb');
fd.storbinary('STOR /Public/sergey/spectroPlots/' + spectroPlotName2_png,fi);
fi.close();

fi = open('spectroPlotName.txt','rb');
fd.storbinary('STOR /Public/sergey/spectroPlots/spectroPlotName.txt',fi);
fi.close();

fd.quit();
