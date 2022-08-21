# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 01:30:35 2020

@author: Arduser
"""

import sys
from PyQt5 import QtWidgets, QtCore, QtGui
import NCtrlGUI
from NCtlLib import *
#from pysnmp.hlapi.asyncore import *
#from pysnmp.proto.api import v2c


class SRH36AntCtrl(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        NCtrlGUI.initUI(self)
        #self.timer_1sec = self.startTimer(1000, QtCore.Qt.VeryCoarseTimer)
    
    def closeEvent(self, event):
        reply = QtWidgets.QMessageBox.question(self, 'Confirm Message',
            "Are you sure to quit?", QtWidgets.QMessageBox.Yes |
            QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
        if reply == QtWidgets.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()
    
    def setMode(self,IP,cmd):
        errorIndication, errorStatus, errorIndex, varBinds = next(
                        setCmd(SnmpEngine(),
                       CommunityData('public', mpModel=0),
                       UdpTransportTarget((IP, 161), retries = 0, timeout = 0.5),
                       ContextData(),
                       ObjectType(ObjectIdentity(OID_BSU+OID_setMode),Integer(cmd))
                       )
                )
        #Ant = list(dAntIP.keys())[list(dAntIP.values()).index(IP)] # получаем номер АП по IP
        if errorIndication:
            self.log.moveCursor(QtGui.QTextCursor.Start)
            self.log.setTextColor(QtGui.QColor("red"))# Вывод в поле log текст красным
            self.log.insertPlainText(format(IP)+ ' - ' + format(errorIndication)+ '\n')
            print(IP, '-' , errorIndication)
            #print(errorIndication)
        elif errorStatus:
            print(IP, '-', '%s at %s' % (errorStatus.prettyPrint(),errorIndex and varBinds[int(errorIndex) - 1][0] or '?'))
            #print('%s at %s' % (errorStatus.prettyPrint(),errorIndex and varBinds[int(errorIndex) - 1][0] or '?'))
        else:
            for varBind in varBinds:
                self.log.moveCursor(QtGui.QTextCursor.Start)
                self.log.setTextColor(QtGui.QColor("black"))
                self.log.insertPlainText(format(IP)+' - '+format(list(dSetAntMode.keys())[list(dSetAntMode.values()).index(varBind[1])])+'\n')
                print(IP, '-', list(dSetAntMode.keys())[list(dSetAntMode.values()).index(varBind[1])])


#    def cbFun(self, snmpEngine, sendRequestHandle, errorIndication,
#              errorStatus, errorIndex, varBindTable, cbCtx):
#        authData, transportTarget = cbCtx
#        fromIP = transportTarget.transportAddr[0]
#   #    print('%s via %s' % (authData, transportTarget))
#        if errorIndication:
#            print(fromIP, errorIndication)
#            return
#       elif errorStatus:
#           print('%s at %s' % (errorStatus.prettyPrint(),
#                                errorIndex and varBindTable[-1][int(errorIndex) - 1][0] or '?'))
#           return
#       else:
##                #print (fromIP, varBindRow[0], '=' , varBindRow[1])
 #              for varBind in varBindRow:
 #                  print(fromIP, ' = '.join([x.prettyPrint() for x in varBind]))
 #          return True  # continue table retrieval

 #   def setAsyncCmd(self, targets):
 #       snmpEngine = SnmpEngine()
 #       # Submit initial GETNEXT requests and wait for responses
 #       for authData, transportTarget, varBinds in targets:
 #           setCmd(snmpEngine, authData, transportTarget, ContextData(),
 #                   *varBinds, cbFun=self.cbFun, cbCtx=(authData, transportTarget))
 #       
  #      snmpEngine.transportDispatcher.runDispatcher()

    def setAntCmd(self,IP,OID,cmd):
        errorIndication, errorStatus, errorIndex, varBinds = next(
                       setCmd(SnmpEngine(),
                       CommunityData('public', mpModel=0),
                       UdpTransportTarget((IP, 161), retries = 0, timeout = 0.5),
                       ContextData(),
                       ObjectType(ObjectIdentity(OID), cmd)
                       )
                )
        #Ant = list(dAntIP.keys())[list(dAntIP.values()).index(IP)] # получаем номер АП по IP
        if errorIndication:
            self.log.moveCursor(QtGui.QTextCursor.Start)
            self.log.setTextColor(QtGui.QColor("red"))# Вывод в поле log текст красным
            self.log.insertPlainText(format(IP)+ ' - ' + format(errorIndication)+ '\n')
            print(IP, '-' , errorIndication)
        elif errorStatus:
            print(IP, '-','%s at %s' % (errorStatus.prettyPrint(),
                                    errorIndex and varBinds[int(errorIndex) - 1][0] or '?'))
        else:
            for varBind in varBinds:
                #Ant = list(dAntIP.keys())[list(dAntIP.values()).index(IP)] # получаем номер АП по IP
                if (OID == OID_BSU + OID_setMode):
                    self.log.moveCursor(QtGui.QTextCursor.Start)
                    self.log.setTextColor(QtGui.QColor("black"))
                    self.log.insertPlainText(format(IP)+' - '        +format(list(dSetAntMode.keys())[list(dSetAntMode.values()).index(varBind[1])])+'\n')
                    print(IP, '-',                                           list(dSetAntMode.keys())[list(dSetAntMode.values()).index(varBind[1])])
                elif(OID == OID_MM + OID_cmdSetAz):
                    self.log.moveCursor(QtGui.QTextCursor.Start)
                    self.log.setTextColor(QtGui.QColor("black"))
                    self.log.insertPlainText(format(IP)+'- Azimuth '  +format(list(dSetCommand.keys())[list(dSetCommand.values()).index(varBind[1])])+'\n')
                    print(IP,                           '- Azimuth ',         list(dSetCommand.keys())[list(dSetCommand.values()).index(varBind[1])])
                elif(OID == OID_MM + OID_cmdSetEl):
                    self.log.moveCursor(QtGui.QTextCursor.Start)
                    self.log.setTextColor(QtGui.QColor("black"))
                    self.log.insertPlainText(format(IP)+'- Elevation '+format(list(dSetCommand.keys())[list(dSetCommand.values()).index(varBind[1])])+'\n')
                    print(IP,                           '- Elevation ',       list(dSetCommand.keys())[list(dSetCommand.values()).index(varBind[1])])

    def setAntAngle(self,IP,OID,Angle):
        HexAngleValue = '9f7804' + hex(struct.unpack('<I', struct.pack('<f', Angle))[0])[2:] # упаковываем float to Hex для отправки через snmp 
        errorIndication, errorStatus, errorIndex, varBinds = next(
                        setCmd(SnmpEngine(),
                       CommunityData('public', mpModel=0),
                       UdpTransportTarget((IP, 161), retries = 0, timeout = 0.5),
                       ContextData(),
                       ObjectType(ObjectIdentity(OID),Opaque(hexValue = HexAngleValue))
                       )
                )
        #Ant = list(dAntIP.keys())[list(dAntIP.values()).index(IP)] # получаем номер АП по IP
        if errorIndication:
            self.log.moveCursor(QtGui.QTextCursor.Start)
            self.log.setTextColor(QtGui.QColor("red"))# Вывод в поле log текст красным
            self.log.insertPlainText(format(IP)+ ' - ' + format(errorIndication)+ '\n')
            print(IP, '-',errorIndication)
        elif errorStatus:
            print(IP, '-','%s at %s' % (errorStatus.prettyPrint(),
                                    errorIndex and varBinds[int(errorIndex) - 1][0] or '?'))
        else:
            for varBind in varBinds:
                value = str(varBind)[len(str(varBind))-8:len(str(varBind))-0]# выборка подобрана экспериментально
                angle = (round(struct.unpack('!f', bytes.fromhex(value))[0],4))
                if (OID == OID_MM + OID_setAngleAz):
                    self.log.moveCursor(QtGui.QTextCursor.Start)
                    self.log.setTextColor(QtGui.QColor("black"))
                    self.log.insertPlainText(format(IP)+' write Azimuth = '  + format(angle)+'\n')
                    print(IP, 'установить значение Азимута = ', angle)
                if (OID == OID_MM + OID_setAngleEl):
                    self.log.moveCursor(QtGui.QTextCursor.Start)
                    self.log.setTextColor(QtGui.QColor("black"))
                    self.log.insertPlainText(format(IP)+' write Elevation = '+ format(angle)+'\n')
                    print(IP, 'установить значение Угла места = ', angle)
#    def onClickbtnAntNorth(self): # Обработчик сигналов нажатия со всех кнопок антенн Севера
#        sender = self.sender()# Ловим сигнал и определяем имя обьекта источника сигнала
#       btnSenderName = sender.objectName()
#        
#    def onClickbtnAntWest(self): # Обработчик сигналов нажатия со всех кнопок антенн Запада
#       sender = self.sender()# Ловим сигнал и определяем имя обьекта источника сигнала
#       btnSenderName = sender.objectName()
#   def onClickbtnAntEast(self): # Обработчик сигналов нажатия со всех кнопок антенн Запада
#        sender = self.sender()# Ловим сигнал и определяем имя обьекта источника сигнала
#        btnSenderName = sender.objectName()
#    def onClickbtnAntFarEast(self): # Обработчик сигналов нажатия со всех кнопок антенн Запада
#        sender = self.sender()# Ловим сигнал и определяем имя обьекта источника сигнала
#        btnSenderName = sender.objectName()
#            
#    def onClickbtnName(self): # Обработчик сигналов нажатия со всех кнопок антенн
#        sender = self.sender()# Ловим сигнал и определяем имя обьекта источника сигнала
#       btnSenderName = sender.objectName()
#        if (self.lst_btn_name[int(btnSenderName)].isChecked() == False):
#                self.btnSRH.setChecked(False)
#        if (int(btnSenderName)<9): 
#            self.lst_btn_direction[0].setChecked(False)
#            if (self.lst_btn_name[0].isChecked()&self.lst_btn_name[1].isChecked()&self.lst_btn_name[2].isChecked()&self.lst_btn_name[3].isChecked()&self.lst_btn_name[4].isChecked()&self.lst_btn_name[5].isChecked()&self.lst_btn_name[6].isChecked()&self.lst_btn_name[7].isChecked()&self.lst_btn_name[8].isChecked()):
#                self.lst_btn_direction[0].setChecked(True)
#                if (self.lst_btn_direction[0].isChecked()&self.lst_btn_direction[1].isChecked()&self.lst_btn_direction[2].isChecked()&self.lst_btn_direction[3].isChecked()):# выделить север при всех нажатых группах
#                    self.btnSRH.setChecked(True)
#        if (9<=int(btnSenderName)<17): 
#            self.lst_btn_direction[1].setChecked(False)
#            if (self.lst_btn_name[0+9].isChecked()&self.lst_btn_name[1+9].isChecked()&self.lst_btn_name[2+9].isChecked()&self.lst_btn_name[3+9].isChecked()&self.lst_btn_name[4+9].isChecked()&self.lst_btn_name[5+9].isChecked()&self.lst_btn_name[6+9].isChecked()&self.lst_btn_name[7+9].isChecked()):
#                self.lst_btn_direction[1].setChecked(True)
#                if (self.lst_btn_direction[0].isChecked()&self.lst_btn_direction[1].isChecked()&self.lst_btn_direction[2].isChecked()&self.lst_btn_direction[3].isChecked()):# выделить север при всех нажатых группах
#                    self.btnSRH.setChecked(True)
#        if (17<=int(btnSenderName)<25): 
#            self.lst_btn_direction[2].setChecked(False)
#            if (self.lst_btn_name[0+9+8].isChecked()&self.lst_btn_name[1+9+8].isChecked()&self.lst_btn_name[2+9+8].isChecked()&self.lst_btn_name[3+9+8].isChecked()&self.lst_btn_name[4+9+8].isChecked()&self.lst_btn_name[5+9+8].isChecked()&self.lst_btn_name[6+9+8].isChecked()&self.lst_btn_name[7+9+8].isChecked()):
#                self.lst_btn_direction[2].setChecked(True)
#               if (self.lst_btn_direction[0].isChecked()&self.lst_btn_direction[1].isChecked()&self.lst_btn_direction[2].isChecked()&self.lst_btn_direction[3].isChecked()):# выделить север при всех нажатых группах
#                   self.btnSRH.setChecked(True)
#        if (25<=int(btnSenderName)<33): 
#            self.lst_btn_direction[3].setChecked(False)
#            if (self.lst_btn_name[0+9+8+8].isChecked()&self.lst_btn_name[1+9+8+8].isChecked()&self.lst_btn_name[2+9+8+8].isChecked()&self.lst_btn_name[3+9+8+8].isChecked()&self.lst_btn_name[4+9+8+8].isChecked()&self.lst_btn_name[5+9+8+8].isChecked()&self.lst_btn_name[6+9+8+8].isChecked()&self.lst_btn_name[7+9+8+8].isChecked()):
#                self.lst_btn_direction[3].setChecked(True)
#                if (self.lst_btn_direction[0].isChecked()&self.lst_btn_direction[1].isChecked()&self.lst_btn_direction[2].isChecked()&self.lst_btn_direction[3].isChecked()):# выделить север при всех нажатых группах
#                    self.btnSRH.setChecked(True)
    def onClickbtnSRH(self): # Выбор антенн и плечей
        #self.setAsyncCmd()# удалить!!! только для тестов
        if (self.btnSRH.isChecked()):
            for i in range(len(self.lst_btn_ant_West)):
                self.lst_btn_ant_West[i].setChecked(True) # снимем выделение антенн Запада при входе в режим Север
                self.lst_btn_ant_West[i].setVisible(False) # убираем кнопки антенн Запада
            for i in range(len(self.lst_btn_ant_East)): 
                self.lst_btn_ant_East[i].setChecked(True) # снимем выделение антенн Востока при выходе из режима управления
                self.lst_btn_ant_East[i].setVisible(False) # убираем кнопки антенн Востока
            for i in range(len(self.lst_btn_ant_FarEast)): 
                self.lst_btn_ant_FarEast[i].setChecked(True) # снимем выделение антенн Востока при выходе из режима управления
                self.lst_btn_ant_FarEast[i].setVisible(False) # убираем кнопки антенн Востока
            for i in range(len(self.lst_btn_ant_North)): 
                self.lst_btn_ant_North[i].setChecked(True) # снимем выделение антенн Севера при выходе из режима управления
                self.lst_btn_ant_North[i].setVisible(False) # убираем кнопки антенн Севера               
            #for i in range(len(self.lst_btn_name)): self.lst_btn_name[i].setChecked(True)
            for i in range(len(self.lst_btn_direction)): self.lst_btn_direction[i].setChecked(True)
        else: 
            #for i in range(len(self.lst_btn_name)): self.lst_btn_name[i].setChecked(False)
            for i in range(len(self.lst_btn_direction  )): self.lst_btn_direction[i].setChecked(False)
            for i in range(len(self.lst_btn_ant_North  )): self.lst_btn_ant_North[i].setChecked(False) #
            for i in range(len(self.lst_btn_ant_West   )): self.lst_btn_ant_West[i].setChecked(False) #
            for i in range(len(self.lst_btn_ant_East   )): self.lst_btn_ant_East[i].setChecked(False) #
            for i in range(len(self.lst_btn_ant_FarEast)): self.lst_btn_ant_FarEast[i].setChecked(False) # 
    def onClickbtnNorth(self):
        for i in range(3): self.lst_btn_direction[i+1].setChecked(False) # убираем выделение с остальных плечей
        #self.lst_btn_direction[0].setChecked()
        if (self.lst_btn_direction[0].isChecked()):
            for i in range(len(self.lst_btn_ant_West)):
                self.lst_btn_ant_West[i].setChecked(False) # снимем выделение антенн Запада при входе в режим Север
                self.lst_btn_ant_West[i].setVisible(False) # убираем кнопки антенн Запада
            for i in range(len(self.lst_btn_ant_East)): 
                self.lst_btn_ant_East[i].setChecked(False) # снимем выделение антенн Востока при выходе из режима управления
                self.lst_btn_ant_East[i].setVisible(False) # убираем кнопки антенн Востока
            for i in range(len(self.lst_btn_ant_FarEast)): 
                self.lst_btn_ant_FarEast[i].setChecked(False) # снимем выделение антенн Востока при выходе из режима управления
                self.lst_btn_ant_FarEast[i].setVisible(False) # убираем кнопки антенн Востока
            for i in range(len(self.lst_btn_ant_North)): 
                self.lst_btn_ant_North[i].setVisible(True) # показываем кнопки антенн Севера
                self.lst_btn_ant_North[i].setChecked(True) #
        else:
            for i in range(len(self.lst_btn_ant_North)): 
                self.lst_btn_ant_North[i].setChecked(False) # снимем выделение антенн Севера при выходе из режима управления
                #self.lst_btn_ant_North[i].setVisible(False) # убираем кнопки антенн Севера
        self.btnSRH.setChecked(False) # снимаем выделение СРГ 3-6
    def onClickbtnWest(self):
        for i in (0,2,3): self.lst_btn_direction[i].setChecked(False) # убираем выделение с остальных плечей
        if (self.lst_btn_direction[1].isChecked()):
            for i in range(len(self.lst_btn_ant_North)): 
                self.lst_btn_ant_North[i].setChecked(False) # снимем выделение антенн Севера при выходе из режима управления
                self.lst_btn_ant_North[i].setVisible(False) # убираем кнопки антенн Севера
            for i in range(len(self.lst_btn_ant_East)): 
                self.lst_btn_ant_East[i].setChecked(False) # снимем выделение антенн Востока при выходе из режима управления
                self.lst_btn_ant_East[i].setVisible(False) # убираем кнопки антенн Востока
            for i in range(len(self.lst_btn_ant_FarEast)): 
                self.lst_btn_ant_FarEast[i].setChecked(False) # снимем выделение антенн Востока при выходе из режима управления
                self.lst_btn_ant_FarEast[i].setVisible(False) # убираем кнопки антенн Востока
            for i in range(len(self.lst_btn_ant_West)): 
                self.lst_btn_ant_West[i].setVisible(True) # показываем кнопки антенн Севера
                self.lst_btn_ant_West[i].setChecked(True) #
        else:
            for i in range(len(self.lst_btn_ant_West)): 
                self.lst_btn_ant_West[i].setChecked(False) # снимем выделение антенн Запада при выходе из режима управления
                #self.lst_btn_ant_West[i].setVisible(False) # убираем кнопки антенн Запада
        self.btnSRH.setChecked(False) # снимаем выделение СРГ 3-6
    
    def onClickbtnEast(self):
        for i in (0,1,3): self.lst_btn_direction[i].setChecked(False) # убираем выделение с остальных плечей
        if (self.lst_btn_direction[2].isChecked()):
            for i in range(len(self.lst_btn_ant_North)): 
                self.lst_btn_ant_North[i].setChecked(False) # снимем выделение антенн Севера при выходе из режима управления
                self.lst_btn_ant_North[i].setVisible(False) # убираем кнопки антенн Севера
            for i in range(len(self.lst_btn_ant_West)): 
                self.lst_btn_ant_West[i].setChecked(False) # снимем выделение антенн Запада при выходе из режима управления
                self.lst_btn_ant_West[i].setVisible(False) # убираем кнопки антенн Запада
            for i in range(len(self.lst_btn_ant_FarEast)): 
                self.lst_btn_ant_FarEast[i].setChecked(False) # снимем выделение антенн Востока при выходе из режима управления
                self.lst_btn_ant_FarEast[i].setVisible(False) # убираем кнопки антенн Востока
            for i in range(len(self.lst_btn_ant_East)): 
                self.lst_btn_ant_East[i].setVisible(True) # показываем кнопки антенн Севера
                self.lst_btn_ant_East[i].setChecked(True) #
        else:
            for i in range(len(self.lst_btn_ant_East)): 
                self.lst_btn_ant_East[i].setChecked(False) # снимем выделение антенн Востока при выходе из режима управления
                #self.lst_btn_ant_East[i].setVisible(False) # убираем кнопки антенн Востока
        self.btnSRH.setChecked(False) # снимаем выделение СРГ 3-6
    
    def onClickbtnFarEast(self):
        for i in (0,1,2): self.lst_btn_direction[i].setChecked(False) # убираем выделение с остальных плечей
        if (self.lst_btn_direction[3].isChecked()):
            for i in range(len(self.lst_btn_ant_North)): 
                self.lst_btn_ant_North[i].setChecked(False) # снимем выделение антенн Севера при выходе из режима управления
                self.lst_btn_ant_North[i].setVisible(False) # убираем кнопки антенн Севера
            for i in range(len(self.lst_btn_ant_West)): 
                self.lst_btn_ant_West[i].setChecked(False) # снимем выделение антенн Запада при выходе из режима управления
                self.lst_btn_ant_West[i].setVisible(False) # убираем кнопки антенн Запада
            for i in range(len(self.lst_btn_ant_East)): 
                self.lst_btn_ant_East[i].setChecked(False) # снимем выделение антенн Востока при выходе из режима управления
                self.lst_btn_ant_East[i].setVisible(False) # убираем кнопки антенн Востока
            for i in range(len(self.lst_btn_ant_FarEast)): 
                self.lst_btn_ant_FarEast[i].setVisible(True) # показываем кнопки антенн Севера
                self.lst_btn_ant_FarEast[i].setChecked(True) #
        else:
            for i in range(len(self.lst_btn_ant_FarEast)): 
                self.lst_btn_ant_FarEast[i].setChecked(False) # снимем выделение антенн Востока при выходе из режима управления
                #self.lst_btn_ant_FarEast[i].setVisible(False) # убираем кнопки антенн Востока
        self.btnSRH.setChecked(False) # снимаем выделение СРГ 3-6

    def onClickbtnSunTracking(self):
        self.btnAzStop.setEnabled(False)# Disable buttons due to needs to stop SUN tracking
        self.btnAzLeft.setEnabled(False)
        self.btnAzRight.setEnabled(False)
        self.btnElStop.setEnabled(False)
        self.btnElUp.setEnabled(False)
        self.btnElDown.setEnabled(False)
        self.spinAzAngle.setDisabled(True);
#        targetIP=[]
#        targets = []
        for antNorth in dAntNorthIP: 
            if (self.lst_btn_ant_North[list(dAntNorthIP.keys()).index(antNorth)].isChecked()): #если антенна выбрана
                self.setMode(dAntNorthIP.get(antNorth),dSetAntMode.get('МСЭ-2'))# берем ее IP и посылыем туда команду
#                targetIP.append(dAntNorthIP.get(antNorth))
        for antWest in dAntWestIP: 
            if (self.lst_btn_ant_West[list(dAntWestIP.keys()).index(antWest)].isChecked()): #если антенна выбрана
                self.setMode(dAntWestIP.get(antWest),dSetAntMode.get('МСЭ-2'))# берем ее IP и посылыем туда команду
#                targetIP.append(dAntWestIP.get(antWest))
        for antEast in dAntEastIP: 
            if (self.lst_btn_ant_East[list(dAntEastIP.keys()).index(antEast)].isChecked()): #если антенна выбрана
                self.setMode(dAntEastIP.get(antEast),dSetAntMode.get('МСЭ-2'))# берем ее IP и посылыем туда команду
#                targetIP.append(dAntEastIP.get(antEast))
        for antFarEast in dAntFarEastIP: 
            if (self.lst_btn_ant_FarEast[list(dAntFarEastIP.keys()).index(antFarEast)].isChecked()): #если антенна выбрана
                self.setMode(dAntFarEastIP.get(antFarEast),dSetAntMode.get('МСЭ-2'))# берем ее IP и посылыем туда команду
#                targetIP.append(dAntFarEastIP.get(antFarEast))
        

#   for IP in targetIP:
#            targets.append(
#            (CommunityData('public'),
#             UdpTransportTarget((IP, 161), retries = 0, timeout = 0.5), 
#            (ObjectType(ObjectIdentity(OID_BSU+OID_setMode),v2c.OctetString('33')),)),)
            #!!!timeout minimal step is 0.5, retries = 1 means send two times, 0 send once
             #by default timeout = 1 second, retries = 5 times
#        targets = tuple(targets)  
        #self.setAsyncCmd(targets)

    def onClickbtnStopMoving(self):
        if (self.btnSunTracking.isChecked):
            self.btnSunTracking.setChecked(False)
        self.btnAzStop.setEnabled(True)
        self.btnAzLeft.setEnabled(True)
        self.btnAzRight.setEnabled(True)
        self.btnElStop.setEnabled(True)
        self.btnElUp.setEnabled(True)
        self.btnElDown.setEnabled(True)
        self.spinAzAngle.setDisabled(False);    
#        for antName in dAntIP: 
#            if (self.lst_btn_name[list(dAntIP.keys()).index(antName)].isChecked()): #если антенна выбрана
#                self.setMode(dAntIP.get(antName),dSetAntMode.get('БСУ выключить'))# берем ее IP и посылыем туда команду
        for antNorth in dAntNorthIP: 
            if (self.lst_btn_ant_North[list(dAntNorthIP.keys()).index(antNorth)].isChecked()): #если антенна выбрана
                self.setMode(dAntNorthIP.get(antNorth),dSetAntMode.get('БСУ выключить'))# берем ее IP и посылыем туда команду
#                targetIP.append(dAntNorthIP.get(antNorth))
        for antWest in dAntWestIP: 
            if (self.lst_btn_ant_West[list(dAntWestIP.keys()).index(antWest)].isChecked()): #если антенна выбрана
                self.setMode(dAntWestIP.get(antWest),dSetAntMode.get('БСУ выключить'))# берем ее IP и посылыем туда команду
#                targetIP.append(dAntWestIP.get(antWest))
        for antEast in dAntEastIP: 
            if (self.lst_btn_ant_East[list(dAntEastIP.keys()).index(antEast)].isChecked()): #если антенна выбрана
                self.setMode(dAntEastIP.get(antEast),dSetAntMode.get('БСУ выключить'))# берем ее IP и посылыем туда команду
#                targetIP.append(dAntEastIP.get(antEast))
        for antFarEast in dAntFarEastIP: 
            if (self.lst_btn_ant_FarEast[list(dAntFarEastIP.keys()).index(antFarEast)].isChecked()): #если антенна выбрана
                self.setMode(dAntFarEastIP.get(antFarEast),dSetAntMode.get('БСУ выключить'))# берем ее IP и посылыем туда команду
#                targetIP.append(dAntFarEastIP.get(antFarEast))
                
    def onClickbtnServicePoint(self):
        if (self.btnSunTracking.isChecked):
            self.btnSunTracking.setChecked(False)
        self.btnAzStop.setEnabled(False)
        self.btnAzLeft.setEnabled(False)
        self.btnAzRight.setEnabled(False)
        self.btnElStop.setEnabled(False)
        self.btnElUp.setEnabled(False)
        self.btnElDown.setEnabled(False)
#        for antName in dAntIP: 
#            if (self.lst_btn_name[list(dAntIP.keys()).index(antName)].isChecked()): #если антенна выбрана
#                self.setAntCmd  (dAntIP.get(antName), OID_BSU + OID_setMode, Integer(dSetAntMode.get('БСУ выключить')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
#                self.setAntAngle(dAntIP.get(antName), OID_MM  + OID_setAngleAz, 0.1234) #Задаем угол по Азимуту
#                self.setAntAngle(dAntIP.get(antName), OID_MM  + OID_setAngleEl, 0.6789) #Задаем угол по Углуместа
#                self.setAntCmd  (dAntIP.get(antName), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
#                self.setAntCmd  (dAntIP.get(antName), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                #print(dAntIP.get(antName))

        for antNorth in dAntNorthIP: 
            if (self.lst_btn_ant_North[list(dAntNorthIP.keys()).index(antNorth)].isChecked()): #если антенна выбрана
                self.setAntCmd  (dAntNorthIP.get(antNorth), OID_BSU + OID_setMode, Integer(dSetAntMode.get('БСУ выключить')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                self.setAntAngle(dAntNorthIP.get(antNorth), OID_MM  + OID_setAngleAz, 0.1234) #Задаем угол по Азимуту
                self.setAntAngle(dAntNorthIP.get(antNorth), OID_MM  + OID_setAngleEl, 0.6789) #Задаем угол по Углуместа
                self.setAntCmd  (dAntNorthIP.get(antNorth), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                self.setAntCmd  (dAntNorthIP.get(antNorth), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                #targetIP.append(dAntNorthIP.get(antNorth))
        for antWest in dAntWestIP: 
            if (self.lst_btn_ant_West[list(dAntWestIP.keys()).index(antWest)].isChecked()): #если антенна выбрана
                self.setAntCmd  (dAntWestIP.get(antWest), OID_BSU + OID_setMode, Integer(dSetAntMode.get('БСУ выключить')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                self.setAntAngle(dAntWestIP.get(antWest), OID_MM  + OID_setAngleAz, 0.1234) #Задаем угол по Азимуту
                self.setAntAngle(dAntWestIP.get(antWest), OID_MM  + OID_setAngleEl, 0.6789) #Задаем угол по Углуместа
                self.setAntCmd  (dAntWestIP.get(antWest), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                self.setAntCmd  (dAntWestIP.get(antWest), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                #targetIP.append(dAntNorthIP.get(antNorth))
        for antEast in dAntEastIP: 
            if (self.lst_btn_ant_East[list(dAntEastIP.keys()).index(antEast)].isChecked()): #если антенна выбрана
                self.setAntCmd  (dAntEastIP.get(antEast), OID_BSU + OID_setMode, Integer(dSetAntMode.get('БСУ выключить')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                self.setAntAngle(dAntEastIP.get(antEast), OID_MM  + OID_setAngleAz, 0.1234) #Задаем угол по Азимуту
                self.setAntAngle(dAntEastIP.get(antEast), OID_MM  + OID_setAngleEl, 0.6789) #Задаем угол по Углуместа
                self.setAntCmd  (dAntEastIP.get(antEast), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                self.setAntCmd  (dAntEastIP.get(antEast), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                #targetIP.append(dAntNorthIP.get(antNorth))
        for antFarEast in dAntFarEastIP: 
            if (self.lst_btn_ant_FarEast[list(dAntFarEastIP.keys()).index(antFarEast)].isChecked()): #если антенна выбрана
                self.setAntCmd  (dAntFarEastIP.get(antFarEast), OID_BSU + OID_setMode, Integer(dSetAntMode.get('БСУ выключить')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                self.setAntAngle(dAntFarEastIP.get(antFarEast), OID_MM  + OID_setAngleAz, 0.1234) #Задаем угол по Азимуту
                self.setAntAngle(dAntFarEastIP.get(antFarEast), OID_MM  + OID_setAngleEl, 0.6789) #Задаем угол по Углуместа
                self.setAntCmd  (dAntFarEastIP.get(antFarEast), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                self.setAntCmd  (dAntFarEastIP.get(antFarEast), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                #targetIP.append(dAntNorthIP.get(antNorth))
    def shiftPointing(self, antIP, Axis, delta):
        Angles = getAntAngle(antIP)# Angles list: [0]- Azimuth, [1] - Elevation
        if (Axis == 0):# 0 - Azimuth, 1 - Elevation
            OID_setAngle = OID_setAngleAz
            OID_cmdSet = OID_cmdSetAz
        elif (Axis == 1):
            OID_setAngle = OID_setAngleEl
            OID_cmdSet = OID_cmdSetEl
        else: print("Axes has a wrong value at shiftPointing()")
        if type(Angles)==list: # если получен угол(есть ответ), устанавливаем значение и отправляем ОПУ на заданный угол
            self.setAntAngle(antIP, OID_MM  + OID_setAngle, Angles[Axis] + delta) #Задаем угол по Азимуту
            self.setAntCmd  (antIP, OID_MM  + OID_cmdSet, Integer(dSetCommand.get('Set Angle')))#
    
    def onClickAzForw(self):
        delta = shiftDelta
        Axis = 0# 0-Azimuth
        for antNorth in dAntNorthIP: 
            if (self.lst_btn_ant_North[list(dAntNorthIP.keys()).index(antNorth)].isChecked()): #если антенна выбрана
                antIP = dAntNorthIP.get(antNorth)
                self.shiftPointing(antIP, Axis, delta)
        for antWest in dAntWestIP: 
            if (self.lst_btn_ant_West[list(dAntWestIP.keys()).index(antWest)].isChecked()): #если антенна выбрана
                antIP = dAntWestIP.get(antWest)
                self.shiftPointing(antIP, Axis, delta)
        for antEast in dAntEastIP: 
            if (self.lst_btn_ant_East[list(dAntEastIP.keys()).index(antEast)].isChecked()): #если антенна выбрана
                antIP = dAntEastIP.get(antEast)
                self.shiftPointing(antIP, Axis, delta)
        for antFarEast in dAntFarEastIP: 
            if (self.lst_btn_ant_FarEast[list(dAntFarEastIP.keys()).index(antFarEast)].isChecked()): #если антенна выбрана
                antIP = dAntFarEastIP.get(antFarEast)
                self.shiftPointing(antIP, Axis, delta)

    def onClickAzBack(self):
        delta = - shiftDelta
        Axis = 0# 0-Elevation, 1-Azimuth
        for antNorth in dAntNorthIP: 
            if (self.lst_btn_ant_North[list(dAntNorthIP.keys()).index(antNorth)].isChecked()): #если антенна выбрана
                antIP = dAntNorthIP.get(antNorth)
                self.shiftPointing(antIP,Axis, delta)
        for antWest in dAntWestIP: 
            if (self.lst_btn_ant_West[list(dAntWestIP.keys()).index(antWest)].isChecked()): #если антенна выбрана
                antIP = dAntWestIP.get(antWest)
                self.shiftPointing(antIP,Axis, delta)
        for antEast in dAntEastIP: 
            if (self.lst_btn_ant_East[list(dAntEastIP.keys()).index(antEast)].isChecked()): #если антенна выбрана
                antIP = dAntEastIP.get(antEast)
                self.shiftPointing(antIP,Axis, delta)
        for antFarEast in dAntFarEastIP: 
            if (self.lst_btn_ant_FarEast[list(dAntFarEastIP.keys()).index(antFarEast)].isChecked()): #если антенна выбрана
                antIP = dAntFarEastIP.get(antFarEast)
                self.shiftPointing(antIP,Axis, delta)

    def onClickAzStop(self):
        for antNorth in dAntNorthIP:
            if (self.lst_btn_ant_North[list(dAntNorthIP.keys()).index(antNorth)].isChecked()): #если антенна выбрана
                self.setAntCmd(dAntNorthIP.get(antNorth), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Waiting')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
        for antWest in dAntWestIP: 
            if (self.lst_btn_ant_West[list(dAntWestIP.keys()).index(antWest)].isChecked()): #если антенна выбрана
                self.setAntCmd(dAntWestIP.get(antWest), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Waiting')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
        for antEast in dAntEastIP: 
            if (self.lst_btn_ant_East[list(dAntEastIP.keys()).index(antEast)].isChecked()): #если антенна выбрана
                self.setAntCmd(dAntEastIP.get(antEast), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Waiting')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
        for antFarEast in dAntFarEastIP:
            if (self.lst_btn_ant_FarEast[list(dAntFarEastIP.keys()).index(antFarEast)].isChecked()): #если антенна выбрана
                self.setAntCmd(dAntFarEastIP.get(antFarEast), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Waiting')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду

    def onClickElUp(self):
        delta = shiftDelta
        Axis = 1# 0-Elevation, 1-Azimuth
        for antNorth in dAntNorthIP: 
            if (self.lst_btn_ant_North[list(dAntNorthIP.keys()).index(antNorth)].isChecked()): #если антенна выбрана
                antIP = dAntNorthIP.get(antNorth)
                self.shiftPointing(antIP,Axis, delta)
        for antWest in dAntWestIP: 
            if (self.lst_btn_ant_West[list(dAntWestIP.keys()).index(antWest)].isChecked()): #если антенна выбрана
                antIP = dAntWestIP.get(antWest)
                self.shiftPointing(antIP,Axis, delta)
        for antEast in dAntEastIP: 
            if (self.lst_btn_ant_East[list(dAntEastIP.keys()).index(antEast)].isChecked()): #если антенна выбрана
                antIP = dAntEastIP.get(antEast)
                self.shiftPointing(antIP,Axis, delta)
        for antFarEast in dAntFarEastIP: 
            if (self.lst_btn_ant_FarEast[list(dAntFarEastIP.keys()).index(antFarEast)].isChecked()): #если антенна выбрана
                antIP = dAntFarEastIP.get(antFarEast)
                self.shiftPointing(antIP,Axis, delta)

    def onClickElDown(self):
        delta = - shiftDelta
        Axis = 1# 0-Elevation, 1-Azimuth
        for antNorth in dAntNorthIP: 
            if (self.lst_btn_ant_North[list(dAntNorthIP.keys()).index(antNorth)].isChecked()): #если антенна выбрана
                antIP = dAntNorthIP.get(antNorth)
                self.shiftPointing(antIP,Axis, delta)
        for antWest in dAntWestIP: 
            if (self.lst_btn_ant_West[list(dAntWestIP.keys()).index(antWest)].isChecked()): #если антенна выбрана
                antIP = dAntWestIP.get(antWest)
                self.shiftPointing(antIP,Axis, delta)
        for antEast in dAntEastIP: 
            if (self.lst_btn_ant_East[list(dAntEastIP.keys()).index(antEast)].isChecked()): #если антенна выбрана
                antIP = dAntEastIP.get(antEast)
                self.shiftPointing(antIP,Axis, delta)
        for antFarEast in dAntFarEastIP: 
            if (self.lst_btn_ant_FarEast[list(dAntFarEastIP.keys()).index(antFarEast)].isChecked()): #если антенна выбрана
                antIP = dAntFarEastIP.get(antFarEast)
                self.shiftPointing(antIP,Axis, delta)

    def onClickElStop(self):
        for antNorth in dAntNorthIP:
            if (self.lst_btn_ant_North[list(dAntNorthIP.keys()).index(antNorth)].isChecked()): #если антенна выбрана
                self.setAntCmd(dAntNorthIP.get(antNorth), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Waiting')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
        for antWest in dAntWestIP: 
            if (self.lst_btn_ant_West[list(dAntWestIP.keys()).index(antWest)].isChecked()): #если антенна выбрана
                self.setAntCmd(dAntWestIP.get(antWest), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Waiting')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
        for antEast in dAntEastIP: 
            if (self.lst_btn_ant_East[list(dAntEastIP.keys()).index(antEast)].isChecked()): #если антенна выбрана
                self.setAntCmd(dAntEastIP.get(antEast), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Waiting')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
        for antFarEast in dAntFarEastIP:
            if (self.lst_btn_ant_FarEast[list(dAntFarEastIP.keys()).index(antFarEast)].isChecked()): #если антенна выбрана
                self.setAntCmd(dAntFarEastIP.get(antFarEast), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Waiting')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
    
    def onClickAzSet(self):
        for antNorth in dAntNorthIP:
            if (self.lst_btn_ant_North[list(dAntNorthIP.keys()).index(antNorth)].isChecked()): #если антенна выбрана
                self.setAntAngle(dAntNorthIP.get(antNorth), OID_MM  + OID_setAngleAz, self.spinAzAngle.value()+0.00000001) #Задаем угол по Азимуту
                self.setAntCmd  (dAntNorthIP.get(antNorth), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
        for antWest in dAntWestIP: 
            if (self.lst_btn_ant_West[list(dAntWestIP.keys()).index(antWest)].isChecked()): #если антенна выбрана
                self.setAntAngle(dAntWestIP.get(antWest), OID_MM  + OID_setAngleAz, self.spinAzAngle.value()+0.00000001) #Задаем угол по Азимуту
                self.setAntCmd  (dAntWestIP.get(antWest), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
        for antEast in dAntEastIP: 
            if (self.lst_btn_ant_East[list(dAntEastIP.keys()).index(antEast)].isChecked()): #если антенна выбрана
                self.setAntAngle(dAntEastIP.get(antEast), OID_MM  + OID_setAngleAz, self.spinAzAngle.value()+0.00000001) #Задаем угол по Азимуту
                self.setAntCmd  (dAntEastIP.get(antEast), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
        for antFarEast in dAntFarEastIP:
            if (self.lst_btn_ant_FarEast[list(dAntFarEastIP.keys()).index(antFarEast)].isChecked()): #если антенна выбрана
                self.setAntAngle(dAntFarEastIP.get(antFarEast), OID_MM  + OID_setAngleAz, self.spinAzAngle.value()+0.00000001) #Задаем угол по Азимуту
                self.setAntCmd  (dAntFarEastIP.get(antFarEast), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
   
    def onClickElSet(self):
        for antNorth in dAntNorthIP:
            if (self.lst_btn_ant_North[list(dAntNorthIP.keys()).index(antNorth)].isChecked()): #если антенна выбрана
                self.setAntAngle(dAntNorthIP.get(antNorth), OID_MM  + OID_setAngleEl, self.spinElAngle.value()+0.00000001) #Задаем угол по Азимуту
                self.setAntCmd  (dAntNorthIP.get(antNorth), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
        for antWest in dAntWestIP: 
            if (self.lst_btn_ant_West[list(dAntWestIP.keys()).index(antWest)].isChecked()): #если антенна выбрана
                self.setAntAngle(dAntWestIP.get(antWest), OID_MM  + OID_setAngleEl, self.spinElAngle.value()+0.00000001) #Задаем угол по Азимуту
                self.setAntCmd  (dAntWestIP.get(antWest), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
        for antEast in dAntEastIP: 
            if (self.lst_btn_ant_East[list(dAntEastIP.keys()).index(antEast)].isChecked()): #если антенна выбрана
                self.setAntAngle(dAntEastIP.get(antEast), OID_MM  + OID_setAngleEl, self.spinElAngle.value()+0.00000001) #Задаем угол по Азимуту
                self.setAntCmd  (dAntEastIP.get(antEast), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
        for antFarEast in dAntFarEastIP:
            if (self.lst_btn_ant_FarEast[list(dAntFarEastIP.keys()).index(antFarEast)].isChecked()): #если антенна выбрана
                self.setAntAngle(dAntFarEastIP.get(antFarEast), OID_MM  + OID_setAngleEl, self.spinElAngle.value()+0.00000001) #Задаем угол по Азимуту
                self.setAntCmd  (dAntFarEastIP.get(antFarEast), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду

        
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    SRH3 = SRH36AntCtrl()
    SRH3.show()
    #exit(app.exec_())
    sys.exit(app.exec_())
    