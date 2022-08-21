# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 13:39:59 2019
@author: A.Gubin
Script for rebooting/swiching SNR-SMART-DIN by SNR-ERD-2.
IP adresses: 
    10.1.2.50 - 0.erd-2s.rao.istp.ac.ru
    10.1.2.51 - 1.erd-2s.rao.istp.ac.ru
    10.1.2.52 - 2.erd-2s.rao.istp.ac.ru
    10.1.2.53 - 3.erd-2s.rao.istp.ac.ru
    10.1.2.54 - 4.erd-2s.rao.istp.ac.ru
    10.1.2.55 - 5.erd-2s.rao.istp.ac.ru
Used pins: 
    10(+5v) and 8(DO1) for rebooting, default delay 3 sec
    10(+5v) and 9(DO2) for switching, Integer(1) is OFF, Integer(2) is ON
Commands:
    for rebooting - OID 1.3.6.1.4.1.40418.2.2.2.1,  Value Integer(1)
    for switching - OID 1.3.6.1.4.1.40418.2.3.2.1,  Value Integer(1) is OFF
                                                    Value Integer(0) is ON
""" 
#!/usr/bin/python
#from pysnmp.hlapi import *
import pysnmp.hlapi
g = pysnmp.hlapi.setCmd(pysnmp.hlapi.SnmpEngine()
           , pysnmp.hlapi.CommunityData('public', mpModel=1)
           , pysnmp.hlapi.UdpTransportTarget(('0.erd-2s.rao.istp.ac.ru', 161)) # 10.1.2.50
           , pysnmp.hlapi.ContextData()
           , pysnmp.hlapi.ObjectType(pysnmp.hlapi.ObjectIdentity('1.3.6.1.4.1.40418.2.2.2.1'), pysnmp.hlapi.Integer(1)) #1 = new value
           )
errorIndication, errorStatus, errorIndex, varBinds = next(g)

print(errorIndication, varBinds)