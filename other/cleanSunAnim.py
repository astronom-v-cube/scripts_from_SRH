#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 01:20:56 2018

@author: sergey
"""
import os, fnmatch;
from optparse import OptionParser;
import pylab as PL
import matplotlib.image as mpimg
import matplotlib.animation as animation

parser = OptionParser();
parser.add_option("-f","--file", dest="filePath", default='/home/svlesovoi/Pictures/20170424/png');
parser.add_option("-e","--ext", dest="fileExt", default='png');
parser.add_option("-s","--size", dest="figureSize", default='12');
parser.add_option("-o","--output", dest="outputPath", default='/home/svlesovoi/srh_20170424_7500_clean.mp4');
(cl_options,cl_args) = parser.parse_args();
path = cl_options.filePath;
outputPath = cl_options.outputPath;

fig, ax = PL.subplots(nrows=1, figsize=(12,12))
ax.axis('off')

def findFiles(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                if os.path.getsize(os.path.join(root,basename)) > 2880:
                    result.append(os.path.join(root,basename));
    return result;

picNames = findFiles(path,'*.png')
picNames.sort();

def updateSunPicture(picIndex):
    img=mpimg.imread(picNames[picIndex])
    imgplot = ax.imshow(img)
    return imgplot, 


ani = animation.FuncAnimation(fig, updateSunPicture, frames=len(picNames), blit=True, repeat=False, interval=100)
ani.save(outputPath)

PL.show()
