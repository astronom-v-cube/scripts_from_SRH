#!/usr/bin/python
# -*- coding: utf-8 -*-
import os;
import pyfits as  FITS;
import matplotlib.pyplot as PLOT;
import numpy;
from optparse import OptionParser;

parser = OptionParser();
parser.add_option("-f","--file",dest="fileName",default="v252f.fits");
parser.add_option("-s","--sample",dest="sampleNumber",default="0");
(cl_options,cl_args) = parser.parse_args();
sample = cl_options.sampleNumber;
inFile = cl_options.fileName;

os.chdir('C:/Users/Sergey/py');
hdulist = FITS.open(inFile);
vis_matrix = hdulist[8].data.field('FLUX');
vis = vis_matrix[sample,:].reshape(2,4,64,2);

PLOT.plot(vis[0,0,:,0]);
PLOT.plot(vis[1,0,:,0]);
PLOT.plot(numpy.sqrt(numpy.power(vis[0,0,:,0],2.) + numpy.power(vis[1,0,:,0],2.)));
PLOT.show();
