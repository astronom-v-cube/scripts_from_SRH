#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 09:07:47 2020

@author: svlesovoi
"""
import numpy as NP

class CSoldatFile():
    def __init__(self, soldatFilePath):
        soldatFileHandler = open(soldatFilePath, encoding='ISO-8859-1')
        soldat_data_text = soldatFileHandler.readlines()
        soldat_data_text = soldat_data_text[10234:]

        self.soldat_decl = NP.zeros((len(soldat_data_text)))
        self.soldat_dDdt = NP.zeros((len(soldat_data_text)))
        self.soldat_d2Ddt2 = NP.zeros((len(soldat_data_text)))
        self.soldat_noon = NP.zeros((len(soldat_data_text)))
        self.soldat_dHdt = NP.zeros((len(soldat_data_text)))
        self.soldat_d2Hdt2 = NP.zeros((len(soldat_data_text)))
        self.soldat_radius = NP.zeros((len(soldat_data_text)))
        self.soldat_pAngle = NP.zeros((len(soldat_data_text)))
        self.soldat_bAngle = NP.zeros((len(soldat_data_text)))
        self.soldat_day2index = {}

        for i in range(len(soldat_data_text)):
            data_strs = soldat_data_text[i].split(chr(0xB3))
            decl_text = data_strs[4].split(':')
            noon_text = data_strs[5].split(':')
            if (decl_text[0][0] != '-'):
                self.soldat_decl[i] = float(decl_text[0]) + float(decl_text[1])/60. + float(decl_text[2])/3600.
            else:
                self.soldat_decl[i] = float(decl_text[0]) - float(decl_text[1])/60. - float(decl_text[2])/3600.
            self.soldat_noon[i] = float(noon_text[0])*3600 + float(noon_text[1])*60. + float(noon_text[2])
            self.soldat_radius[i] = float(data_strs[6])
            self.soldat_pAngle[i] = float(data_strs[7])
            self.soldat_bAngle[i] = float(data_strs[8])
            self.soldat_dDdt[i] = NP.deg2rad(float(data_strs[10]) / 3600.) / 3600.
            self.soldat_day2index[data_strs[1] + '-' + ('%02d'%int(data_strs[2])) + '-' + ('%02d'%int(data_strs[3]))] = i

        self.soldat_decl = NP.deg2rad(self.soldat_decl)
        self.soldat_d2Ddt2 = (self.soldat_dDdt[1:] - self.soldat_dDdt[0:-1]) / (24*3600.)
        self.soldat_hourAngle = NP.deg2rad(self.soldat_noon*15)
        self.soldat_dHdt = 2.*NP.pi / (24.*3600 + (self.soldat_noon[1:] - self.soldat_noon[0:-1]))
        self.soldat_d2Hdt2 = (self.soldat_dHdt[1:] - self.soldat_dHdt[0:-1]) / (24*3600.)
