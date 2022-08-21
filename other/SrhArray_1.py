#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 04:14:14 2018

@author: sergey
"""
import numpy as NP

class SrhArray():
    def __init__(self, parent = 0):
        self.adjacentSpacing = 4.9 # meters
        self.phi = 0.903338787600965 # SRH latitude
        self.EastFirst = 65
        self.EastLast = 128
        self.WestFirst = 1
        self.WestLast = 64
        
    def baseline2uvw(self, hourAngle, declination, antennaA, antennaB):
        if ((antennaA >= 1 and antennaA <= 128) and (antennaB >= 129 and antennaB <= 256)):
            base = NP.array([192.5 - antennaB, antennaA - 64.5,0.])
        elif ((antennaB >= 1 and antennaB <= 128) and (antennaA >= 129 and antennaA <= 256)):
            base = NP.array([192.5 - antennaA ,antennaB - 64.5,0.])
        elif ((antennaA >= 1 and antennaA <= 128) and (antennaB >= 1 and antennaB <= 128)):
            base = NP.array([0.,antennaA - antennaB,0.])
        elif ((antennaA >= 129 and antennaA <= 256) and (antennaB >= 129 and antennaB <= 256)):
            base = NP.array([antennaA - antennaB,0.,0.])
        
        base *= self.adjacentSpacing;
        
        phi_operator = NP.array([
            [-NP.sin(self.phi), 0., NP.cos(self.phi)],
            [0., 1., 0.],
            [NP.cos(self.phi), 0., NP.sin(self.phi)]
            ])
    
        uvw_operator = NP.array([
            [ NP.sin(hourAngle),		 NP.cos(hourAngle),		0.	  ],
            [-NP.sin(declination) * NP.cos(hourAngle),   NP.sin(declination) * NP.sin(hourAngle), NP.cos(declination)], 
            [ NP.cos(declination) * NP.cos(hourAngle), -NP.cos(declination) * NP.sin(hourAngle), NP.sin(declination)]  
            ])
    
        return NP.dot(uvw_operator, NP.dot(phi_operator, base))
        