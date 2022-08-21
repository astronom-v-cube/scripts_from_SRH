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
        self.dishDiameter = 1.8 # meters
        self.phi = 0.903338787600965 # SRH latitude
        self.EastFirst = 65
        self.EastLast = 128
        self.WestFirst = 1
        self.WestLast = 64
        self.SouthFirst = 129
        self.SouthLast = 192
        self.NorthFirst = 193
        self.NorthLast = 256
        self.phi_operator = NP.array([
            [-NP.sin(self.phi), 0., NP.cos(self.phi)],
            [0., 1., 0.],
            [NP.cos(self.phi), 0., NP.sin(self.phi)]
            ])
        
    def baseline2uvw(self, hourAngle, declination, antennaA, antennaB):
        if ((antennaA >= self.WestFirst and antennaA <= self.EastLast) and (antennaB >= self.SouthFirst and antennaB <= self.NorthLast)):
            base = NP.array([self.SouthLast + .5 - antennaB, antennaA - self.WestLast - .5,0.])
        elif ((antennaB >= self.WestFirst and antennaB <= self.EastLast) and (antennaA >= self.SouthFirst and antennaA <= self.NorthLast)):
            base = NP.array([self.SouthLast + .5 - antennaA ,antennaB - self.WestLast - .5,0.])
        elif ((antennaA >= self.WestFirst and antennaA <= self.EastLast) and (antennaB >= self.WestFirst and antennaB <= self.EastLast)):
            base = NP.array([0.,antennaA - antennaB,0.])
        elif ((antennaA >= self.SouthFirst and antennaA <= self.NorthLast) and (antennaB >= self.SouthFirst and antennaB <= self.NorthLast)):
            base = NP.array([antennaA - antennaB,0.,0.])
        
        base *= self.adjacentSpacing;
        
        uvw_operator = NP.array([
            [ NP.sin(hourAngle),		 NP.cos(hourAngle),		0.	  ],
            [-NP.sin(declination) * NP.cos(hourAngle),   NP.sin(declination) * NP.sin(hourAngle), NP.cos(declination)], 
            [ NP.cos(declination) * NP.cos(hourAngle), -NP.cos(declination) * NP.sin(hourAngle), NP.sin(declination)]  
            ])
    
        return NP.dot(uvw_operator, NP.dot(self.phi_operator, base))
    
    def isShadowed(self, baseline):
        return baseline - self.dishDiameter
    
    def correctBaselineForShadowing(self, baseline):
        return 0.5*(self.dishDiameter - baseline)
    
        