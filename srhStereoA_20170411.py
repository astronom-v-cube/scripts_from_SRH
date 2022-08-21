# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 01:47:01 2016

@author: Sergey
"""
import pylab as PL;
from optparse import OptionParser;
from astropy.io import fits
from scipy.io import readsav
import urllib.request
import numpy as NP

def hhmm_format(t, pos):
  hh = int(t / 3600.)
  t -= hh*3600.
  mm = int(t / 60.)
  t -= mm*60.
  ss = int(t)
  return '%02d:%02d:%02d' % (hh,mm,ss);

def hm_format(t, pos):
  hh = int(t / 60.)
  t -= hh*60.
  mm = int(t)
  return '%02d:%02d' % (hh,mm);

def freq_format(f, pos):
    if (f < 368):
        return '%3.1f' % (stereoA['frequencies'][int(f)] * 1e-3)
    else:
        return ''

parser = OptionParser();
parser.add_option("-d", "--date", dest="inDate",   default="20170411");
(cl_options,cl_args) = parser.parse_args();
requestedDate = cl_options.inDate;

srhCpName = 'srh_cp_' + requestedDate + '.fits'
#srhURL = 'http://archive.badary.iszf.irk.ru/SRH/' + srhCpName
srhURL = 'http://archive.rao.istp.ac.ru/SRH/' + srhCpName
urllib.request.urlretrieve(srhURL, srhCpName)

stereoName = 'swaves_average_' + requestedDate + '_a.sav'
stereoURL = 'https://solar-radio.gsfc.nasa.gov/data/stereo/new_summary/' + requestedDate[0:4] + '/' + stereoName
urllib.request.urlretrieve(stereoURL, stereoName)

cpFits = fits.open(srhCpName)
stereoA = readsav(stereoName)

cpTime = cpFits[2].data['time']
cpI = cpFits[2].data['I']
cpV = cpFits[2].data['V']
freqs = cpTime.shape[0]

mcpI_time = cpI.mean(0)
mcpI_freq = cpI.mean(1)

goesFits = fits.open('go1320170411.fits')
goesTime = goesFits[2].data['time']
goesFlux = goesFits[2].data['Flux']

fig = PL.figure(figsize=(20,16))
fig.suptitle(cpFits[0].header['DATE-OBS'])
pl0 = PL.subplot(211)
pl0.xaxis.set_major_formatter(PL.FuncFormatter(hhmm_format))
pl0.xaxis.set_major_locator(PL.MultipleLocator(3600))
#pl0.set_ylim(-0.01,0.15)
pl0.set_ylim(-0.01,NP.round(1.2*cpI.max(),2))
pl0.set_xlim(0,10*3600)

for i in range(freqs):
    pl0.plot(cpTime[i], cpI[i], color=PL.cm.gray(int(i/freqs*200)),label='%d MHz'%cpFits[1].data['frequencies'][i])
    pl0.plot(cpTime[i], cpV[i], color=PL.cm.gray(int(i/freqs*200)))
pl0.plot(goesTime[0,:],goesFlux[0,:,0]*3e4,color='red')
pl0.plot(goesTime[0,:],goesFlux[0,:,1]*3e4,color='blue')
pl0.legend()

pl0.grid()
stereoTime = NP.linspace(0,1440*60,1440)
pl0.plot(stereoTime, stereoA['spectrum'].T[100]*1e-3,'.', color='red')
pl0.plot(stereoTime, stereoA['spectrum'].T[200]*1e-3,'.', color='yellow')
pl0.plot(stereoTime, stereoA['spectrum'].T[300]*1e-3,'.', color='green')
pl0.plot(stereoTime, stereoA['spectrum'].T[100]*1e-3, linewidth=0.5, color=PL.cm.gray(50))
pl0.plot(stereoTime, stereoA['spectrum'].T[200]*1e-3, linewidth=0.5, color=PL.cm.gray(100))
pl0.plot(stereoTime, stereoA['spectrum'].T[300]*1e-3, linewidth=0.5, color=PL.cm.gray(150))

pl1 = PL.subplot(212)
pl1.xaxis.set_major_formatter(PL.FuncFormatter(hm_format))
pl1.xaxis.set_major_locator(PL.MultipleLocator(60))
pl1.yaxis.set_major_formatter(PL.FuncFormatter(freq_format))
pl1.set_xlim(0,600)
pl1.set_ylim(1,350)

pl1.imshow(stereoA['spectrum'].T,origin='lower',vmin=0, vmax=10., aspect=1/3)
#pl1.plot(cpTime[0]/60, 3000*mcpI_time, color='red', linewidth=0.5)
pl1.plot([cpTime[0,0]/60,cpTime[0,-1]/60], [100,100], color='red', linewidth=1.)
pl1.plot([cpTime[0,0]/60,cpTime[0,-1]/60], [200,200], color='yellow', linewidth=1.)
pl1.plot([cpTime[0,0]/60,cpTime[0,-1]/60], [300,300], color='green', linewidth=1.)
#pl1.set_yscale('log')
pl1.set_xlabel('UT')
pl1.set_ylabel('MHz')

pl1.plot(goesTime[0,:]/60,goesFlux[0,:,0]*1e8,color='white')
pl1.plot(goesTime[0,:]/60,goesFlux[0,:,1]*1e8,color='white')
goesFits.close()
