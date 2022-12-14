# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 03:56:57 2020

@author: Arduser
"""
from pysnmp.hlapi import *
from pysnmp import proto
import struct
import NCtrlGUI

OID_BSU = '1.3.6.1.4.1.19707.12.1.2.1.'
OID_MM  = '1.3.6.1.4.1.19707.12.2.2.1.'
OID_PDU = '1.3.6.1.4.1.19707.12.4.2.1.'
OID_MP  = '1.3.6.1.4.1.19707.12.5.2.1.'
#Controling
OID_setMode ="4.9.0" #unit display="Команда" format="integer" name="WorkCmd"
dSetAntMode = {'БСУ выключить': 0, 'МСЭ-1': 17, 'МСЭ-2':33, 'BROADCAST-1':49, 'BROADCAST-2':65, 'внешнее управление':161} # словарь заполнения режима работы ОПУ
OID_setAngleAz   = '4.2.1.0'#unit ="Азимут заданный"     format="float" unitname="°" precision="4" min="-180.5" max="180.5" name="param4021"
OID_setAngleEl   = '4.2.2.0'#unit ="Угол места заданный" format="float" unitname="°" precision="4" min="-0.5" max="90.5" name="param4022"
OID_cmdSetAz     = '4.1.1.0'#unit ="Команда, азимут"     format="integer" name="param4011"
OID_cmdSetEl     = '4.1.2.0'#unit ="Команда, угол места" format="integer" name="param4012"
OID_getAngleAz   = '3.2.1.0'#param ="Азимут"     format="float" unitname="°" precision="4" name="param3021"
OID_getAngleEl   = '3.2.2.0'#param ="Угол места" format="float" unitname="°" precision="4" name="param3022"

shiftDelta = 1.5 #Shift antenna point to this value 

dSetCommand={'Waiting'           :0,
            'Set Angle'          :17,
            'Move by table'      :49,
            'Angle stabilization':65,
            'Left move'          :257,
            'Right move'         :273,
            'Encoder calibration':289}

#dAntIP = {'C00' : '10.1.1.200',
#          'N01' : '10.1.13.1' ,'N02':'10.1.13.2' ,'N03':'10.1.13.3' ,'N04':'10.1.13.4' ,'N05':'10.1.13.5',
#          'N06' : '10.1.13.6' ,'N07':'10.1.13.7' ,'N08':'10.1.13.8' ,'N09':'10.1.13.9' ,'N10':'10.1.13.10',
#          'N11' : '10.1.13.11','N12':'10.1.13.12','N13':'10.1.13.13','N14':'10.1.13.14','N15':'10.1.13.15',
#          'N16' : '10.1.13.16','N17':'10.1.13.17','N18':'10.1.13.18','N19':'10.1.13.19','N20':'10.1.13.20',
#          'N21' : '10.1.13.21','N22':'10.1.13.22','N23':'10.1.13.23','N24':'10.1.13.24','N25':'10.1.13.25',
#          'N26' : '10.1.13.26','N27':'10.1.13.27','N28':'10.1.13.28','N29':'10.1.13.29','N30':'10.1.13.30'} #,'N31' : '10.1.13.31'}

dAntNorthIP = {'C00' : '10.1.13.100',
          'N01' : '10.1.13.31' ,'N02':'10.1.13.32' ,'N03':'10.1.13.33' ,'N04':'10.1.13.34' ,'N05':'10.1.13.35',
          'N06' : '10.1.13.36' ,'N07':'10.1.13.37' ,'N08':'10.1.13.38' ,'N09':'10.1.13.39' ,'N10':'10.1.13.40',
          'N11' : '10.1.13.41' ,'N12':'10.1.13.42' ,'N13':'10.1.13.43' ,'N14':'10.1.13.44' ,'N15':'10.1.13.45',
          'N16' : '10.1.13.46' ,'N17':'10.1.13.47' ,'N18':'10.1.13.48' ,'N19':'10.1.13.49' ,'N20':'10.1.13.50',
          'N21' : '10.1.13.51' ,'N22':'10.1.13.52' ,'N23':'10.1.13.53' ,'N24':'10.1.13.54' ,'N25':'10.1.13.55',
          'N26' : '10.1.13.56' ,'N27':'10.1.13.57' ,'N28':'10.1.13.58' ,'N29':'10.1.13.59' ,'N30':'10.1.13.60',
          'N31' : '10.1.13.61' ,'N32':'10.1.13.62'} #
dAntWestIP = {
          'W01' : '10.1.13.201','W02':'10.1.13.202','W03':'10.1.13.203','W04':'10.1.13.204','W05':'10.1.13.205',
          'W06' : '10.1.13.206','W07':'10.1.13.207','W08':'10.1.13.208','W09':'10.1.13.209','W10':'10.1.13.210',
          'W11' : '10.1.13.211','W12':'10.1.13.212','W13':'10.1.13.213','W14':'10.1.13.214','W15':'10.1.13.215',
          'W16' : '10.1.13.216','W17':'10.1.13.217','W18':'10.1.13.218','W19':'10.1.13.219','W20':'10.1.13.220',
          'W21' : '10.1.13.221','W22':'10.1.13.222','W23':'10.1.13.223','W24':'10.1.13.224','W25':'10.1.13.225',
          'W26' : '10.1.13.226','W27':'10.1.13.227','W28':'10.1.13.228','W29':'10.1.13.229','W30':'10.1.13.230',
          'W31' : '10.1.13.231','W32':'10.1.13.232'} #
dAntEastIP = {
          'E01' : '10.1.13.101','E02':'10.1.13.102','E03':'10.1.13.103','E04':'10.1.13.104','E05':'10.1.13.105',
          'E06' : '10.1.13.106','E07':'10.1.13.107','E08':'10.1.13.108','E09':'10.1.13.109','E10':'10.1.13.110',
          'E11' : '10.1.13.111','E12':'10.1.13.112','E13':'10.1.13.113','E14':'10.1.13.114','E15':'10.1.13.115',
          'E16' : '10.1.13.116','E17':'10.1.13.117','E18':'10.1.13.118','E19':'10.1.13.119','E20':'10.1.13.120',
          'E21' : '10.1.13.121','E22':'10.1.13.122','E23':'10.1.13.123','E24':'10.1.13.124','E25':'10.1.13.125',
          'E26' : '10.1.13.126','E27':'10.1.13.127','E28':'10.1.13.128','E29':'10.1.13.129','E30':'10.1.13.130',
          'E31' : '10.1.13.131','E32':'10.1.13.132'} #
dAntFarEastIP = {
          'E33' : '10.1.13.133','E34':'10.1.13.134','E35':'10.1.13.135','E36':'10.1.13.136','E37':'10.1.13.137',
          'E38' : '10.1.13.138','E39':'10.1.13.139','E40':'10.1.13.140','E41':'10.1.13.141','E42':'10.1.13.142',
          'E43' : '10.1.13.143','E44':'10.1.13.144','E45':'10.1.13.145','E46':'10.1.13.146','E47':'10.1.13.147',
          'E48' : '10.1.13.148','E49':'10.1.13.149','E50':'10.1.13.150','E51':'10.1.13.151','E52':'10.1.13.152',
          'E53' : '10.1.13.153','E54':'10.1.13.154','E55':'10.1.13.155','E56':'10.1.13.156','E57':'10.1.13.157',
          'E58' : '10.1.13.158','E59':'10.1.13.159','E60':'10.1.13.160','E61':'10.1.13.161','E62':'10.1.13.162',
          'E63' : '10.1.13.163','E64':'10.1.13.164'} #
                
def getAntAngle(IP):# получение действительных углов из АП
    Angles = [] # список получаемых углов
    errorIndication, errorStatus, errorIndex, varBinds = next(
            getCmd(SnmpEngine(),
                   CommunityData('public', mpModel=0),
                   UdpTransportTarget((IP, 161), retries = 0, timeout = 0.5),
                   ContextData(),
                   ObjectType(ObjectIdentity(OID_MM+OID_getAngleAz)),
                   ObjectType(ObjectIdentity(OID_MM+OID_getAngleEl)),
                   lookupMib=False,
                   lexicographicMode=False)
            )
    if errorIndication:
        print(IP, '-',errorIndication)
    elif errorStatus:
        print(IP, '-','%s at %s' % (errorStatus.prettyPrint(), errorIndex and varBinds[int(errorIndex) - 1][0] or '?'))
    else:
        for varBind in varBinds:
            value = str(varBind)[len(str(varBind))-11:len(str(varBind))-3] # вытаскиваем из строки нужное значение 4 байта of HEX
            Angles.append((round(struct.unpack('!f', bytes.fromhex(value))[0],4)))
        return Angles