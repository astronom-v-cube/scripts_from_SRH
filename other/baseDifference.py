#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 02:56:07 2022

@author: sergey_lesovoi
"""

import numpy as NP
import pylab as PL
from scipy.stats import linregress
from scipy import signal
import scipy.optimize as opt

def fitPhi(h, C, dx, dy):
    return C + dx*NP.cos(h) - dy*NP.sin(h)

def fitDPhi(h, dx, dy):
    return -dx*NP.sin(h) - dy*NP.cos(h)

length = phaseEdit1224.srhFits.visLcp.shape[1]
scanAccum = 120
scanStep = 10
slopeNumber = (length - scanAccum) // scanStep
freqNumber = 10
freqs = phaseEdit.srhFits.freqList * 1e3

lcpPhaseSlope = NP.zeros((slopeNumber, freqNumber,47))
rcpPhaseSlope = NP.zeros((slopeNumber, freqNumber,47))
for ss in range(slopeNumber):
    for ff in range(freqNumber):
        for bb in range(47):
            phaT = NP.unwrap(NP.angle(phaseEdit1224.srhFits.visLcp[ff,ss*scanStep:ss*scanStep + scanAccum,11775 + bb]))
            phaSlope, interc, A,B,C = linregress(NP.linspace(0,scanAccum-1,scanAccum),phaT)
            lcpPhaseSlope[ss,ff,bb] = phaSlope
            phaT = NP.unwrap(NP.angle(phaseEdit1224.srhFits.visRcp[ff,ss*scanStep:ss*scanStep + scanAccum,11775 + bb]))
            phaSlope, interc, A,B,C = linregress(NP.linspace(0,scanAccum-1,scanAccum),phaT)
            rcpPhaseSlope[ss,ff,bb] = phaSlope
        
refBase = 9
lcpBaseRel = NP.zeros((slopeNumber,47))
rcpBaseRel = NP.zeros((slopeNumber,47))
for ss in range(slopeNumber):
    for bb in range(47):
        lcpBaseRel[ss,bb] = (lcpPhaseSlope[ss,:,bb]/lcpPhaseSlope[ss,:,refBase]).mean()
        rcpBaseRel[ss,bb] = (rcpPhaseSlope[ss,:,bb]/rcpPhaseSlope[ss,:,refBase]).mean()
    
PL.figure()
PL.title(phaseEdit1224.srhFits.dateObs)
#PL.ylim(.95,1.05)
PL.ylabel('relative length')
PL.xlabel('base')
PL.grid()
PL.plot(0.5*(rcpBaseRel.T + lcpBaseRel.T),'.')

baseLengthLcp = NP.zeros((length,47))
baseLengthRcp = NP.zeros((length,47))
for bb in range(47):
    for ff in range(freqNumber):
        baseLengthLcp[:,bb] += NP.unwrap(NP.angle(phaseEdit1224.srhFits.visLcp[ff,:,11775 + bb]))/(2*NP.pi*freqs[ff])*3e8
        baseLengthRcp[:,bb] += NP.unwrap(NP.angle(phaseEdit1224.srhFits.visRcp[ff,:,11775 + bb]))/(2*NP.pi*freqs[ff])*3e8
    baseLengthLcp[:,bb] /= freqNumber
    baseLengthRcp[:,bb] /= freqNumber

baseLengthI = 0.5*(baseLengthLcp + baseLengthRcp)

hAngle = []
for ss in range(length):
    hAngle.append(phaseEdit1224.srhFits.getHourAngle(ss))
dH = hAngle[1] - hAngle[2]

baseSlope = []
for bb in range(47):
#    bSlope, interc, A,B,C = linregress(NP.linspace(0,800-1,800),baseLengthI[0:800,bb])
    bSlope, interc, A,B,C = linregress(hAngle,baseLengthI[:,bb])
    baseSlope.append(bSlope)



outDx = []
outDy = []
fitPars = []
fitDPars = []

for base in range(47):
    rBase = baseLengthI[:,base]
    rDBase = NP.gradient(rBase,dH)
    rDDBase = NP.gradient(rDBase,dH)
    
    fitDPar = opt.curve_fit(fitDPhi,hAngle,rDBase[0:800],p0=[0.001,0.001])
    fitPar  = opt.curve_fit(fitPhi, hAngle,rBase[0:800], p0=[0,fitDPar[0][0],fitDPar[0][1]])
    fitPars.append(fitPar[0])
    fitDPars.append(fitDPar[0])
    
    resDXDY = []
    resDDXDDY = []
    Amatr = NP.zeros((2,2))
    for ss in range(800):
        hA = phaseEdit1224.srhFits.getHourAngle(ss)
        Amatr[0,0] =  NP.cos(hA)
        Amatr[0,1] = -NP.sin(hA)
        Amatr[1,0] = -NP.sin(hA)
        Amatr[1,1] = -NP.cos(hA)
        invAmatr = NP.linalg.pinv(Amatr)
        resDXDY.append(invAmatr.dot([rBase[ss],rDBase[ss]]))
    
        Amatr[0,0] = -NP.sin(hA)
        Amatr[0,1] = -NP.cos(hA)
        Amatr[1,0] = -NP.cos(hA)
        Amatr[1,1] =  NP.sin(hA)
        invAmatr = NP.linalg.pinv(Amatr)
        resDDXDDY.append(invAmatr.dot([rDBase[ss],rDDBase[ss]]))
    resDXDY = NP.array(resDXDY)
    resDDXDDY = NP.array(resDDXDDY)
    outDx.append(resDXDY[:,0].mean())
    outDy.append(resDXDY[:,1].mean())

dX = 0.001
dY = 0.005
N = 300
hAng = NP.linspace(.3,.9,N)

dXcos = dX*NP.cos(hAng)
dXsin = dX*NP.sin(hAng)
dYcos = dY*NP.cos(hAng)
dYsin = dY*NP.sin(hAng)

phi = dXcos - dYsin + 0.00 + NP.random.randn(N)*1e-4
# dPhi = -dXsin - dYcos
# ddPhi = -dXcos + dYsin
dPhi = NP.gradient(phi, hAng[1] - hAng[0])
ddPhi = NP.gradient(dPhi, hAng[1] - hAng[0])

resDXDY = []
resDDXDDY = []
Amatr = NP.zeros((2,2))
for ss in range(N):
    Amatr[0,0] =  NP.cos(hAng[ss])
    Amatr[0,1] = -NP.sin(hAng[ss])
    Amatr[1,0] = -NP.sin(hAng[ss])
    Amatr[1,1] = -NP.cos(hAng[ss])
    invAmatr = NP.linalg.pinv(Amatr)
    resDXDY.append(invAmatr.dot([phi[ss],dPhi[ss]]))

    Amatr[0,0] = -NP.sin(hAng[ss])
    Amatr[0,1] = -NP.cos(hAng[ss])
    Amatr[1,0] = -NP.cos(hAng[ss])
    Amatr[1,1] =  NP.sin(hAng[ss])
    invAmatr = NP.linalg.pinv(Amatr)
    resDDXDDY.append(invAmatr.dot([dPhi[ss],ddPhi[ss]]))
resDXDY = NP.array(resDXDY)
resDDXDDY = NP.array(resDDXDDY)


