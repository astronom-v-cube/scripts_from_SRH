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
from optparse import OptionParser

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

def findFits(path, pattern, size):
  result = [];
  for root, dirs, files in os.walk(path):
    for basename in files:
      if fnmatch.fnmatch(basename,pattern):
	if os.path.getsize(os.path.join(root,basename)) >= size:
	  result.append(os.path.join(root,basename));
  return result;

os.environ['DISPLAY'] = ':0.0';

parser = OptionParser();
parser.add_option("-d", "--date", dest="inDate",   default="20120215");
#parser.add_option("-s", "--size", dest="fitsSize", type="int",  default=172800);
parser.add_option("-s", "--size", dest="fitsSize", type="int",  default=2880);
(cl_options,cl_args) = parser.parse_args();
dateName = cl_options.inDate;

fdName = open('aCurrentData.txt','rb');
#dateName = fdName.readline();
fdName.close();

fitName = findFits('/home/serg/SSRT/data/spectro/','*smf*' + dateName + '*.fit',cl_options.fitsSize);
#fitName = findFits('/home/serg/SSRT/data/spectro/201204/','*smf*' + dateName + '*.fit',cl_options.fitsSize);

fitName.sort();

i = 0;
for fName in fitName:
  hdu = FITS.open(fName);

  if i == 0:
    tStartGlobal= DT.datetime.strptime(hdu[0].header['DATE-OBS'] + ' ' + hdu[0].header['TIME-OBS'], '%Y/%m/%d %H:%M:%S.%f');
    freqs	= hdu[1].data.field('FREQUENCY')*1e-9;
    freq_time	= hdu[1].data.field('TIME');
    rcp_data	= hdu[1].data.field('RCP');
    lcp_data	= hdu[1].data.field('LCP');
    freqAmount  = freqs.shape[0];
  else:
    if i == len(fitName) - 1:
      tStopGlobal = DT.datetime.strptime(hdu[0].header['DATE-OBS'] + ' ' + hdu[0].header['TIME-OBS'], '%Y/%m/%d %H:%M:%S.%f');
    freq_time = concatenate(( freq_time ,hdu[1].data.field('TIME')),axis=1);
    rcp_data  = concatenate((rcp_data,hdu[1].data.field('RCP')),axis=1);
    lcp_data  = concatenate((lcp_data,hdu[1].data.field('LCP')),axis=1);
  i = i + 1;
  hdu.close();

i = 0;

L = rcp_data.shape[1];
F = rcp_data.shape[0];

#sky_t0 = 8030;
#sky_t1 = 8100;
sky_t0 = L - 1 - 200;
sky_t1 = L - 1 - 100;
sun_t0 = L - 1 - 4200;
sun_t1 = L - 1 - 4000;
#sun_t0 = 1000;
#sun_t1 = 1500;
#sky_t0 = 11480;
#sky_t1 = 11520;
#sun_t0 = 11590;
#sun_t1 = 11610;

sky_lcp_level = [0.] * F;
sky_rcp_level = [0.] * F;
sun_lcp_level = [0.] * F;
sun_rcp_level = [0.] * F;

for i in range(0,F):
	sky_lcp_level[i] = lcp_data[i,sky_t0:sky_t1].mean();
	sky_rcp_level[i] = rcp_data[i,sky_t0:sky_t1].mean();
	sun_lcp_level[i] = (lcp_data[i,sun_t0:sun_t1] - sky_lcp_level[i]).mean();
	sun_rcp_level[i] = (rcp_data[i,sun_t0:sun_t1] - sky_rcp_level[i]).mean();

fdClbSky = open('sclb_sky_tmp.txt','wb');
fdClbSun = open('sclb_sun_tmp.txt','wb');
for i in range(0,F):
	fdClbSky.write("%.6f %.6f\n" % (sky_lcp_level[i], sky_rcp_level[i]));
	fdClbSun.write("%.6f %.6f\n" % (sun_lcp_level[i], sun_rcp_level[i]));
fdClbSky.close();
fdClbSun.close();

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
fdClbSky = open('sclb_sky_tmp.txt','rb');
for clbSkyLine in fdClbSky.readlines():
  clbSkyRow = clbSkyLine.split();
  for clbSkyWord in clbSkyRow:
    clbSky.append(float(clbSkyWord));
min_sky = array(clbSky).reshape(F,2).transpose();
fdClbSky.close();

clbSun = [];
fdClbSun = open('sclb_sun_tmp.txt','rb');
for clbSunLine in fdClbSun.readlines():
  clbSunRow = clbSunLine.split();
  for clbSunWord in clbSunRow:
    clbSun.append(float(clbSunWord));
mean_sun = array(clbSun).reshape(F,2).transpose();
fdClbSun.close();

for i in range(0,F):
#  lcp_data[i,:] = lcp_data[i,:] - min_sky[0,i];
#  lcp_data[i,:] = lcp_data[i,:]/mean_sun[0,i] * sp_2_24_flux[i];
#  rcp_data[i,:] = rcp_data[i,:] - min_sky[1,i];
#  rcp_data[i,:] = rcp_data[i,:]/mean_sun[1,i] * sp_2_24_flux[i];
  lcp_data[i,:] = lcp_data[i,:] - sky_lcp_level[i];
  lcp_data[i,:] = lcp_data[i,:]/sun_lcp_level[i] * sp_2_24_flux[i];
  rcp_data[i,:] = rcp_data[i,:] - sky_rcp_level[i];
  rcp_data[i,:] = rcp_data[i,:]/sun_rcp_level[i] * sp_2_24_flux[i];

for i in range(1,L):
  if (freq_time[0,i] == -1.):
    freq_time[0,i] = freq_time[0,i - 1];

i_data = (rcp_data + lcp_data)/2.;
v_data = (rcp_data - lcp_data)/2.;

c_list = matplotlib.colors.LinearSegmentedColormap.from_list(LAB.cm.datad['gist_rainbow'], colors=['r','g','b'], N = F);

fig1 = LAB.figure(1,figsize=(20,12));
LAB.subplot(111);

dt_major = 3600;
dt_minor = 600;

sp1 = LAB.subplot(111);
commonTitle = 'Siberian spectropolarimter, ' + tStartGlobal.strftime('%Y %m %d');
for freq in range (0,freqAmount):
  commonTitle += ", %1.2f" % freqs[freq];
commonTitle += " GHz";
LAB.title(commonTitle,fontsize='large');
LAB.ylabel('I, V [s.f.u.]');
LAB.xlabel('UT');
sp1.xaxis.set_major_locator(MultipleLocator(dt_major));
sp1.xaxis.set_major_formatter(LAB.FuncFormatter(hhmm_format));
sp1.xaxis.set_minor_locator(MultipleLocator(dt_minor));
for i in range(0,F):
  if (i != 11 and i != 15):
    LAB.plot(freq_time[0,:],i_data[i,:],color = c_list(float(i)/F));
for i in range(0,F):
  if (i != 11 and i != 15):
    LAB.plot(freq_time[0,:],v_data[i,:],color = c_list(float(i)/F));
LAB.xlim((0.*3600.,10.5*3600.));
LAB.ylim((-100.,2500.));

#sp2 = LAB.subplot(212);
#LAB.ylabel('RCP - LCP [s.f.u.]');
#LAB.xlabel('UT');
#sp2.xaxis.set_major_locator(MultipleLocator(dt_major));
#sp2.xaxis.set_major_formatter(LAB.FuncFormatter(hhmm_format));
#sp2.xaxis.set_minor_locator(MultipleLocator(dt_minor));
#for i in range(0,F):
#  if (i != 11 and i != 15):
#    LAB.plot(freq_time[0,:],v_data[i,:],color = c_list(float(i)/F));
#LAB.xlim((1.25*3600.,2.25*3600.));
#LAB.ylim((-800.,800.));

#LAB.draw();
#LAB.show();
