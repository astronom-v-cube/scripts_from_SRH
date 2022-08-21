#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 08:24:14 2021

@author: svlesovoi
"""
import json
import numpy as NP
import pylab as PL
import os, fnmatch

def findJson(path, pattern):
    result = [];
    for root, dirs, files in os.walk(path):
        for basename in files:
            if fnmatch.fnmatch(basename,pattern):
                result.append(os.path.join(root,basename));
    return result;

telErr = []
jsonFiles = findJson('srh_tracker','*08-12*.json')
jsonFiles.sort()
errAzHist = NP.zeros(250)
errAltHist = NP.zeros(250)
errAzList = {}
errAltList = {}

for jsonFile in jsonFiles:
    telFd = open(jsonFile)
    telLog = json.load(telFd)
    for logRecord in telLog:
        if 'ip' in logRecord:
            telErr.append(logRecord)
            ind = int(telErr[-1]['ip'].split('.')[-1])
            if 'Код ошибки, азимут' in logRecord:
                errAzHist[ind] += 1
                errAzList[telErr[-1]['ip']] = errAzHist[ind]
            elif 'Код ошибки, угол места' in logRecord:
                errAltHist[ind] += 1
                errAltList[telErr[-1]['ip']] = errAltHist[ind]
    telFd.close()

PL.figure()
PL.plot(errAzHist,'.')
#PL.plot([31,31],[0,5000],color='black')
#PL.plot([62,62],[0,5000],color='black')
#PL.plot([101,101],[0,5000],color='black')
#PL.plot([164,164],[0,5000],color='black')
#PL.plot([201,201],[0,5000],color='black')
#PL.plot([232,232],[0,5000],color='black')
#PL.ylim(0,3500)
PL.yscale('log')
PL.xlabel('IP')
#PL.title('SRH 20210717 azimuth errors')
PL.title('SRH azimuth errors')
PL.grid()

PL.figure()
PL.plot(errAltHist,'.')
#PL.plot([31,31],[0,50],color='black')
#PL.plot([62,62],[0,50],color='black')
#PL.plot([101,101],[0,50],color='black')
#PL.plot([164,164],[0,50],color='black')
#PL.plot([201,201],[0,50],color='black')
#PL.plot([232,232],[0,50],color='black')
#PL.ylim(0,35)
PL.xlabel('IP')
#PL.title('SRH 20210717 altitude errors')
PL.title('SRH altitude errors')
PL.yscale('log')
PL.grid()
