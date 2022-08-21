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

def findFits(path, pattern):
  result = [];
  for root, dirs, files in os.walk(path):
    for basename in files:
      if fnmatch.fnmatch(basename,pattern):
	if os.path.getsize(os.path.join(root,basename)) > 2880:
	  result.append(os.path.join(root,basename));
  return result;

def base2uvw(hourAngle, delta, antenna0, antenna1):
  phi = 0.903338787600965;
  
  if ((antenna0 >= 1 and antenna0 <= 128) and (antenna1 >= 129 and antenna1 <= 256)):
    base = array([192.5 - antenna1,antenna0 - 64.5,0.]);
  elif ((antenna1 >= 1 and antenna1 <= 128) and (antenna0 >= 129 and antenna0 <= 256)):
    base = array([192.5 - antenna0,antenna1 - 64.5,0.]);
  elif ((antenna0 >= 1 and antenna0 <= 128) and (antenna1 >= 1 and antenna1 <= 128)):
    base = array([0.,antenna0 - antenna1,0.]);
  elif ((antenna0 >= 129 and antenna0 <= 256) and (antenna1 >= 129 and antenna1 <= 256)):
    base = array([antenna0 - antenna1,0.,0.]);

  base = base*4.9;

  phi_operator = array([
    [-sin(phi),	0.,	cos(phi)],
    [ 0.,	1.,	0.     ],
    [ cos(phi),	0.,	sin(phi)] 
  ]);

  uvw_operator = array([
    [ sin(hourAngle),		 cos(hourAngle),		0.	  ], 
    [-sin(delta)*cos(hourAngle), sin(delta)*sin(hourAngle),	cos(delta)], 
    [ cos(delta)*cos(hourAngle),-cos(delta)*sin(hourAngle),	sin(delta)]  
  ]);

  base_phi = dot(phi_operator,base);
  uvw =  dot(uvw_operator,base_phi);

  return uvw;

os.environ['DISPLAY'] = ':0.0';

parser = OptionParser();
parser.add_option("-d", "--date", dest="inDate",   default="20100522");
parser.add_option("-s", "--size", dest="fitsSize", type="int",  default=2880);
(cl_options,cl_args) = parser.parse_args();
dateName = cl_options.inDate;

fdName = open('/home/serg/SSRT/data/prototype/aCurrentData.txt','rb');
dateName = fdName.readline();
fdName.close();

fitName = findFits('/home/serg/SSRT/data/prototype/','mf_' + dateName + '*.fit')

fitName.sort();

i = 0;
for fName in fitName:
  hdu = FITS.open(fName);
  freqs	= hdu[1].data.field('FREQUENCY');
  sub_freq_time	= hdu[1].data.field('TIME');
  freqAmount  = freqs.shape[0];
  sub_length  = sub_freq_time.shape[1];
  sub_ant_rcp = zeros((10,freqAmount,sub_length),'float');
  sub_ant_lcp = zeros((10,freqAmount,sub_length),'float');
  sub_vis_rcp = zeros((55,freqAmount,sub_length),'complex');
  sub_vis_lcp = zeros((55,freqAmount,sub_length),'complex');

  if i == 0:
    timeOfStart= DT.datetime.strptime(hdu[0].header['DATE-OBS'] + ' ' + hdu[0].header['TIME-OBS'], '%Y/%m/%d %H:%M:%S.%f');
    ant_rcp = zeros((10,freqAmount,sub_length),'float');
    ant_lcp = zeros((10,freqAmount,sub_length),'float');
    vis_rcp = zeros((55,freqAmount,sub_length),'complex');
    vis_lcp = zeros((55,freqAmount,sub_length),'complex');
    freq_time = hdu[1].data.field('TIME');
    for ant in range(0,10):
      ant_rcp[ant,:,:] = hdu[1].data.field(2*ant + 5);
      ant_lcp[ant,:,:] = hdu[1].data.field(2*ant + 6);
    for vis in range(0,55):
      vis_rcp[vis,:,:] = hdu[1].data.field(2*vis + 25);
      vis_lcp[vis,:,:] = hdu[1].data.field(2*vis + 26);
  else:
    freq_time = concatenate(( freq_time ,hdu[1].data.field('TIME')),axis=1);
    for ant in range(0,10):
      sub_ant_rcp[ant,:] = hdu[1].data.field(2*ant + 5);
      sub_ant_lcp[ant,:] = hdu[1].data.field(2*ant + 6);
    ant_rcp = concatenate((ant_rcp,sub_ant_rcp),axis=2);
    ant_lcp = concatenate((ant_lcp,sub_ant_lcp),axis=2);
    for vis in range(0,55):
      sub_vis_rcp[vis,:,:] = hdu[1].data.field(2*vis + 25);
      sub_vis_lcp[vis,:,:] = hdu[1].data.field(2*vis + 26);
    vis_rcp = concatenate((vis_rcp,sub_vis_rcp),axis=2);
    vis_lcp = concatenate((vis_lcp,sub_vis_lcp),axis=2);
  i = i + 1;
  hdu.close();
  sub_ant_rcp = 0;
  sub_ant_lcp = 0;
  sub_vis_rcp = 0;
  sub_vis_lcp = 0;

L = ant_rcp.shape[2];
F = ant_rcp.shape[1];

ant_i = ant_rcp + ant_lcp;
i_data = vis_rcp + vis_lcp;

A1	= 1.;
A2	= 1.;
A3	= 1.;
A126	= 1.;
A127	= 1.;
A128	= 1.;
A129	= 0.;
A130	= 1.;
A131	= 1.;	
A132	= 1.;

An = zeros(10);
An[0]	= A132;
An[1]	= A131;
An[2]	= A1;
An[3]	= A2;
An[4]	= A129;
An[5]	= A130;
An[6]	= A127;
An[7]	= A128;
An[8]	= A3;
An[9]	= A126;

AnAm = outer(An,An);

A_n_lcp 	= zeros((10,F),float);
A_n_rcp 	= zeros((10,F),float);
vis_flux_lcp	= zeros((55,F,L),float);
vis_flux_rcp	= zeros((55,F,L),float);
vis_phas_lcp	= zeros((55,F,L),float);
vis_phas_rcp	= zeros((55,F,L),float);

flux_sun = array([138.,152.,173.,194.,216.]);

fdNorpData = open('norp_data.txt','rb');
norpLines = fdNorpData.readlines();
norpData = [];
for i in norpLines:
	norpData.append([float(j) for j in i.split()]);
norpData = array(norpData);
fdNorpData.close();

norp_freq = norpData[:,0]*1e3;
norp_flux = norpData[:,1];
norp_fit_coef = polyfit(norp_freq,norp_flux,2);

fit_freq = 20e3*arange(1000.)/1000.;
fit_flux = polyval(norp_fit_coef,fit_freq);

flux_sun = zeros(F);
for i in range(0,F):
	flux_sun[i] = fit_flux[where(abs(fit_freq - freqs[i]) < 0.12)];

zero_flux_t = 200;

clbSky = [];
fdClbSky = open('clb_sky_lcp.txt','rb');
for clbSkyLine in fdClbSky.readlines():
  clbSkyRow = clbSkyLine.split();
  for clbSkyWord in clbSkyRow:
    clbSky.append(float(clbSkyWord));
min_sky_lcp = array(clbSky).reshape(F,10).transpose();
fdClbSky.close();

clbSky = [];
fdClbSky = open('clb_sky_rcp.txt','rb');
for clbSkyLine in fdClbSky.readlines():
  clbSkyRow = clbSkyLine.split();
  for clbSkyWord in clbSkyRow:
    clbSky.append(float(clbSkyWord));
min_sky_rcp = array(clbSky).reshape(F,10).transpose();
fdClbSky.close();

clbSun = [];
fdClbSun = open('clb_sun_lcp.txt','rb');
for clbSunLine in fdClbSun.readlines():
  clbSunRow = clbSunLine.split();
  for clbSunWord in clbSunRow:
    clbSun.append(float(clbSunWord));
mean_sun_lcp = array(clbSun).reshape(F,10).transpose();
fdClbSun.close();

clbSun = [];
fdClbSun = open('clb_sun_rcp.txt','rb');
for clbSunLine in fdClbSun.readlines():
  clbSunRow = clbSunLine.split();
  for clbSunWord in clbSunRow:
    clbSun.append(float(clbSunWord));
mean_sun_rcp = array(clbSun).reshape(F,10).transpose();
fdClbSun.close();

for i in range(55):
  for j in range(F):
    vis_flux_lcp[i,j,zero_flux_t:] = abs(vis_lcp[i,j,zero_flux_t:]);
    vis_flux_rcp[i,j,zero_flux_t:] = abs(vis_rcp[i,j,zero_flux_t:]);
    vis_phas_lcp[i,j,zero_flux_t:] = arctan2(imag(vis_lcp[i,j,zero_flux_t:]),real(vis_lcp[i,j,zero_flux_t:]))*180./pi;
    vis_phas_rcp[i,j,zero_flux_t:] = arctan2(imag(vis_rcp[i,j,zero_flux_t:]),real(vis_rcp[i,j,zero_flux_t:]))*180./pi;

for i in range(10):
  for j in range(F):
    vis_flux_lcp[i,j,:]	= vis_flux_lcp[i,j,:] - min_sky_lcp[i,j];
    vis_flux_rcp[i,j,:]	= vis_flux_rcp[i,j,:] - min_sky_rcp[i,j];
    A_n_lcp[i,j]	= flux_sun[j]/mean_sun_lcp[i,j];
    A_n_rcp[i,j]	= flux_sun[j]/mean_sun_rcp[i,j];

B_n_lcp = sqrt(A_n_lcp);
B_n_rcp = sqrt(A_n_lcp);
nm_ = zeros((2,55),int);
ind = 0;
for i in range(10):
  for j in range(i,10):
    if (i == 0):
      nm_[:,ind] = [j,j];
    else:
      nm_[:,ind] = [i-1,j];
    ind = ind + 1;

for freq in range(F):
  B_nm_lcp = outer(B_n_lcp[:,freq],B_n_lcp[:,freq]);
  B_nm_rcp = outer(B_n_rcp[:,freq],B_n_rcp[:,freq]);
  for i in range(55):
    vis_flux_lcp[i,freq,:] = B_nm_lcp[nm_[0,i],nm_[1,i]]*vis_flux_lcp[i,freq,:];
    vis_flux_rcp[i,freq,:] = B_nm_rcp[nm_[0,i],nm_[1,i]]*vis_flux_rcp[i,freq,:];

#-------------------------------------------------------------------------------------------------
lf = transpose(vis_flux_lcp[0,:,:] + vis_flux_rcp[0,:,:])*An[0];
for i in range(1,10):
  lf = lf + transpose(vis_flux_lcp[i,:,:] + vis_flux_rcp[i,:,:])*An[i];
lf = 0.5 * lf / sum(An);

#-------------------------------------------------------------------------------------------------
ind = array([27,38,32,51,49,53,40,22,21]);
mask = array([A1*A2,A2*A3,A1*A3,A126*A127,A127*A128,A126*A128,A129*A130,A130*A131,A129*A131]);

mf = transpose(vis_flux_lcp[ind[0],:,:] + vis_flux_rcp[ind[0],:,:])*mask[0];
for i in range(1,len(ind)):
  mf = mf + transpose(vis_flux_lcp[ind[i],:,:] + vis_flux_rcp[ind[i],:,:])*mask[i];
mf = 0.5 * mf / sum(mask);

#-------------------------------------------------------------------------------------------------
ind = array([42,41,44,43,34,28,46,45,48,47,35,29,24,23,26,25,20,19,16,15,18,17,12,11,33,30,31,39,36,37,54,50,52]);
mask = array([A129*A128, A129*A127, A129*A126, A129*A3, A129*A2, A129*A1, 
	      A130*A128, A130*A127, A130*A126, A130*A3, A130*A2, A130*A1,
	      A131*A128, A131*A127, A131*A126, A131*A3, A131*A2, A131*A1,
	      A132*A128, A132*A127, A132*A126, A132*A3, A132*A2, A132*A1,
	      A1*A126,   A1*A127,   A1*A128,   A2*A126, A2*A127, A2*A128, A3*A126, A3*A127, A3*A128]);

hf = transpose(vis_flux_lcp[ind[0],:,:] + vis_flux_rcp[ind[0],:,:])*mask[0];
for i in range(1,len(ind)):
  hf = hf + transpose(vis_flux_lcp[ind[i],:,:] + vis_flux_rcp[ind[i],:,:])*mask[i];
hf = 0.5 * hf / sum(mask);

dt_major = 3600.;
dt_minor = 900.;
t0_axis_hours = 0.;
t1_axis_hours = 10.;

zero_ind1 = where(abs(freq_time[0,:] - roll(freq_time[0,:],1)) > 10);
zero_ind1 = array(zero_ind1[0]);
zero_ind2 = where(abs(freq_time[0,:] - roll(freq_time[0,:],-1)) > 10);
zero_ind2 = array(zero_ind2[0]);
if (zero_ind1.shape > 0):
  lf[zero_ind1,:] = 0.;
  mf[zero_ind1,:] = 0.;
  hf[zero_ind1,:] = 0.;
  
if (zero_ind2.shape > 0):
  lf[zero_ind2,:] = 0.;
  mf[zero_ind2,:] = 0.;
  hf[zero_ind2,:] = 0.;

c_list = matplotlib.colors.LinearSegmentedColormap.from_list(LAB.cm.datad['gist_rainbow'], colors=['r','g','b'], N = F);

fig1 = LAB.figure(1,figsize=(20,12));
commonTitle = 'Siberian radioheliograph prototype, ' + timeOfStart.strftime('%Y %B %d');
for freq in range (0,freqAmount):
  commonTitle += ", %1.2f" % (freqs[freq]*1e-3);
commonTitle += " GHz";
fig1.suptitle(commonTitle,fontsize='large');

sp1 = fig1.add_subplot(3,1,1);
sp1.set_ylabel('single dish [s.f.u.]');
sp1.set_xlabel('UT');
sp1.xaxis.set_major_locator(MultipleLocator(dt_major));
sp1.xaxis.set_major_formatter(LAB.FuncFormatter(hhmm_format));
sp1.xaxis.set_minor_locator(MultipleLocator(dt_minor));
for i in range(0,F):
  sp1.plot(freq_time[0,:],lf[:,i],color=c_list(float(i)/F));
sp1.set_xlim((t0_axis_hours*3600.,t1_axis_hours*3600.));
sp1.set_ylim((0.,5000.));

sp2 = fig1.add_subplot(3,1,2);
sp2.set_ylabel('low UV correlation [s.f.u.]');
sp2.set_xlabel('UT');
sp2.xaxis.set_major_locator(MultipleLocator(dt_major));
sp2.xaxis.set_major_formatter(LAB.FuncFormatter(hhmm_format));
sp2.xaxis.set_minor_locator(MultipleLocator(dt_minor));
for i in range(0,F):
  sp2.plot(freq_time[0,:],mf[:,i],color=c_list(float(i)/F));
sp2.set_xlim((t0_axis_hours*3600.,t1_axis_hours*3600.));
sp2.set_ylim((0.,2000.));

sp3 = fig1.add_subplot(3,1,3);
sp3.set_ylabel('high UV correlation [s.f.u.]');
sp3.set_xlabel('UT');
sp3.xaxis.set_major_locator(MultipleLocator(dt_major));
sp3.xaxis.set_major_formatter(LAB.FuncFormatter(hhmm_format));
sp3.xaxis.set_minor_locator(MultipleLocator(dt_minor));
for i in range(0,F):
  sp3.plot(freq_time[0,:],hf[:,i],color=c_list(float(i)/F));
sp3.set_xlim((t0_axis_hours*3600.,t1_axis_hours*3600.));
sp3.set_ylim((0.,100.));

corrPlotName = 'fCorrPlot'+ DT.datetime.now().strftime("%Y%m%d");
corrPlotName_png = corrPlotName + '.png';

LAB.savefig(corrPlotName_png);

fdCPlotName = open('corrPlotName.txt','wb');
fdCPlotName.write(corrPlotName_png);
fdCPlotName.close();

fd = ftplib.FTP('192.168.0.1','sergey','jill21');

fi = open(corrPlotName_png,'rb');
fd.storbinary('STOR /Public/sergey/corrPlots/' + corrPlotName_png,fi);
fi.close();

fi = open('corrPlotName.txt','rb');
fd.storbinary('STOR /Public/sergey/corrPlots/corrPlotName.txt',fi);
fi.close();

fd.quit();
