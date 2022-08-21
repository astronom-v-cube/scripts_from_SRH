#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 01:08:57 2020

@author: svlesovoi
"""
import numpy as NP
import pylab as PL

soldat_data_file = open('/home/svlesovoi/SSRT/soldat.txt', encoding='ISO-8859-1')
soldat_data_text = soldat_data_file.readlines()
soldat_data_text = soldat_data_text[10234:]

soldat_decl = NP.zeros((len(soldat_data_text)))
soldat_dDdt = NP.zeros((len(soldat_data_text)))
soldat_noon = NP.zeros((len(soldat_data_text)))
soldat_radius = NP.zeros((len(soldat_data_text)))
soldat_pAngle = NP.zeros((len(soldat_data_text)))
soldat_bAngle = NP.zeros((len(soldat_data_text)))

for i in range(len(soldat_data_text)):
    data_strs = soldat_data_text[i].split(chr(0xB3))
    decl_text = data_strs[4].split(':')
    noon_text = data_strs[5].split(':')
    if (decl_text[0][0] != '-'):
        soldat_decl[i] = float(decl_text[0]) + float(decl_text[1])/60. + float(decl_text[2])/3600.
    else:
        soldat_decl[i] = float(decl_text[0]) - float(decl_text[1])/60. - float(decl_text[2])/3600.
    soldat_noon[i] = float(noon_text[0])*3600 + float(noon_text[1])*60. + float(noon_text[2])
    soldat_radius[i] = float(data_strs[6])
    soldat_pAngle[i] = float(data_strs[7])
    soldat_bAngle[i] = float(data_strs[8])
    soldat_dDdt[i] = NP.deg2rad(float(data_strs[10]) / 3600.) / 3600.

soldat_decl = NP.deg2rad(soldat_decl)
soldat_d2Ddt2 = (soldat_dDdt[1:] - soldat_dDdt[0:-1]) / (24*3600.)
soldat_hourAngle = NP.deg2rad(soldat_noon*15)
soldat_dHdt = 2.*NP.pi / (24.*3600 + (soldat_noon[1:] - soldat_noon[0:-1]))
soldat_d2Hdt2 = (soldat_dHdt[1:] - soldat_dHdt[0:-1]) / (24*3600.)

PL.figure(figsize=(10,20))
PL.suptitle('Declination')
#PL.suptitle('Hour angle')
pl0 = PL.subplot(311)
pl1 = PL.subplot(312)
pl2 = PL.subplot(313)

pl0.set_ylabel('declination [rad]')
#pl1.set_ylabel('dDdt [rad/s]')
pl1.set_ylabel('radius [arcmin]')
#pl2.set_ylabel('[d2Ddt2 rad/s^2]')
pl2.set_ylabel('P angle [deg]')

pl0.plot(soldat_decl)
#pl1.plot(soldat_dDdt)
#pl2.plot(soldat_d2Ddt2)
pl1.plot(soldat_radius)
pl2.plot(soldat_pAngle)

#pl0.set_ylabel('noon [seconds UTC]')
#pl1.set_ylabel('dHdt [rad/s]')
#pl2.set_ylabel('d2Hdt2 [rad/s^2]')
#pl2.set_xlabel('days')
#
#pl0.plot(soldat_noon)
#pl1.plot(soldat_dHdt)
#pl2.plot(soldat_d2Hdt2)

pl0.grid()
pl1.grid()
pl2.grid()
