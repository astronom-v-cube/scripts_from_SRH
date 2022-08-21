# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 01:30:35 2020

@author: Arduser
"""

import sys
from PyQt5 import QtWidgets, QtCore, QtGui
import NCtrlGUI
from NCtlLib import *

class NorthAntCtrl(QtWidgets.QMainWindow):
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
                       UdpTransportTarget((IP, 161)),
                       ContextData(),
                       ObjectType(ObjectIdentity(OID_BSU+OID_setMode),Integer(cmd))
                       )
                )
        Ant = list(dAntIP.keys())[list(dAntIP.values()).index(IP)] # получаем номер АП по IP
        if errorIndication:
            self.log.moveCursor(QtGui.QTextCursor.Start)
            self.log.setTextColor(QtGui.QColor("red"))# Вывод в поле log текст красным
            self.log.insertPlainText(format(Ant)+ ' - ' + format(errorIndication)+ '\n')
            print(Ant, '-' , errorIndication)
            #print(errorIndication)
        elif errorStatus:
            print(Ant, '-', '%s at %s' % (errorStatus.prettyPrint(),errorIndex and varBinds[int(errorIndex) - 1][0] or '?'))
            #print('%s at %s' % (errorStatus.prettyPrint(),errorIndex and varBinds[int(errorIndex) - 1][0] or '?'))
        else:
            for varBind in varBinds:
                self.log.moveCursor(QtGui.QTextCursor.Start)
                self.log.setTextColor(QtGui.QColor("black"))
                self.log.insertPlainText(format(Ant)+' - '+format(list(dSetAntMode.keys())[list(dSetAntMode.values()).index(varBind[1])])+'\n')
                print(Ant, '-', list(dSetAntMode.keys())[list(dSetAntMode.values()).index(varBind[1])])
    
    def setAntCmd(self,IP,OID,cmd):
        errorIndication, errorStatus, errorIndex, varBinds = next(
                       setCmd(SnmpEngine(),
                       CommunityData('public', mpModel=0),
                       UdpTransportTarget((IP, 161)),
                       ContextData(),
                       ObjectType(ObjectIdentity(OID), cmd)
                       )
                )
        Ant = list(dAntIP.keys())[list(dAntIP.values()).index(IP)] # получаем номер АП по IP
        if errorIndication:
            self.log.moveCursor(QtGui.QTextCursor.Start)
            self.log.setTextColor(QtGui.QColor("red"))# Вывод в поле log текст красным
            self.log.insertPlainText(format(Ant)+ ' - ' + format(errorIndication)+ '\n')
            print(Ant, '-' , errorIndication)
        elif errorStatus:
            print(Ant, '-','%s at %s' % (errorStatus.prettyPrint(),
                                    errorIndex and varBinds[int(errorIndex) - 1][0] or '?'))
        else:
            for varBind in varBinds:
                Ant = list(dAntIP.keys())[list(dAntIP.values()).index(IP)] # получаем номер АП по IP
                if (OID == OID_BSU + OID_setMode):
                    self.log.moveCursor(QtGui.QTextCursor.Start)
                    self.log.setTextColor(QtGui.QColor("black"))
                    self.log.insertPlainText(format(Ant)+' - '        +format(list(dSetAntMode.keys())[list(dSetAntMode.values()).index(varBind[1])])+'\n')
                    print(Ant, '-',                                           list(dSetAntMode.keys())[list(dSetAntMode.values()).index(varBind[1])])
                elif(OID == OID_MM + OID_cmdSetAz):
                    self.log.moveCursor(QtGui.QTextCursor.Start)
                    self.log.setTextColor(QtGui.QColor("black"))
                    self.log.insertPlainText(format(Ant)+'- Azimuth '  +format(list(dSetCommand.keys())[list(dSetCommand.values()).index(varBind[1])])+'\n')
                    print(Ant,                           '- Azimuth ',         list(dSetCommand.keys())[list(dSetCommand.values()).index(varBind[1])])
                elif(OID == OID_MM + OID_cmdSetEl):
                    self.log.moveCursor(QtGui.QTextCursor.Start)
                    self.log.setTextColor(QtGui.QColor("black"))
                    self.log.insertPlainText(format(Ant)+'- Elevation '+format(list(dSetCommand.keys())[list(dSetCommand.values()).index(varBind[1])])+'\n')
                    print(Ant,                           '- Elevation ',       list(dSetCommand.keys())[list(dSetCommand.values()).index(varBind[1])])

    def setAntAngle(self,IP,OID,Angle):
        HexAngleValue = '9f7804' + hex(struct.unpack('<I', struct.pack('<f', Angle))[0])[2:] # упаковываем float to Hex для отправки через snmp 
        errorIndication, errorStatus, errorIndex, varBinds = next(
                        setCmd(SnmpEngine(),
                       CommunityData('public', mpModel=0),
                       UdpTransportTarget((IP, 161)),
                       ContextData(),
                       ObjectType(ObjectIdentity(OID),Opaque(hexValue = HexAngleValue))
                       )
                )
        Ant = list(dAntIP.keys())[list(dAntIP.values()).index(IP)] # получаем номер АП по IP
        if errorIndication:
            self.log.moveCursor(QtGui.QTextCursor.Start)
            self.log.setTextColor(QtGui.QColor("red"))# Вывод в поле log текст красным
            self.log.insertPlainText(format(Ant)+ ' - ' + format(errorIndication)+ '\n')
            print(Ant, '-',errorIndication)
        elif errorStatus:
            print(Ant, '-','%s at %s' % (errorStatus.prettyPrint(),
                                    errorIndex and varBinds[int(errorIndex) - 1][0] or '?'))
        else:
            for varBind in varBinds:
                value = str(varBind)[len(str(varBind))-8:len(str(varBind))-0]# выборка подобрана экспериментально
                angle = (round(struct.unpack('!f', bytes.fromhex(value))[0],4))
                if (OID == OID_MM + OID_setAngleAz):
                    self.log.moveCursor(QtGui.QTextCursor.Start)
                    self.log.setTextColor(QtGui.QColor("black"))
                    self.log.insertPlainText(format(Ant)+' write Azimuth = '  + format(angle)+'\n')
                    print(Ant, 'установить значение Азимута = ', angle)
                if (OID == OID_MM + OID_setAngleEl):
                    self.log.moveCursor(QtGui.QTextCursor.Start)
                    self.log.setTextColor(QtGui.QColor("black"))
                    self.log.insertPlainText(format(Ant)+' write Elevation = '+ format(angle)+'\n')
                    print(Ant, 'установить значение Угла места = ', angle)
        
    def onClickbtnName(self): # Обработчик сигналов нажатия со всех кнопок антенн
        sender = self.sender()# Ловим сигнал и определяем имя обьекта источника сигнала
        btnSenderName = sender.objectName()
        if (self.lst_btn_name[int(btnSenderName)].isChecked() == False):
                self.btnNorth.setChecked(False)
        if (int(btnSenderName)<9): 
            self.lst_btn_group[0].setChecked(False)
            if (self.lst_btn_name[0].isChecked()&self.lst_btn_name[1].isChecked()&self.lst_btn_name[2].isChecked()&self.lst_btn_name[3].isChecked()&self.lst_btn_name[4].isChecked()&self.lst_btn_name[5].isChecked()&self.lst_btn_name[6].isChecked()&self.lst_btn_name[7].isChecked()&self.lst_btn_name[8].isChecked()):
                self.lst_btn_group[0].setChecked(True)
                if (self.lst_btn_group[0].isChecked()&self.lst_btn_group[1].isChecked()&self.lst_btn_group[2].isChecked()&self.lst_btn_group[3].isChecked()):# выделить север при всех нажатых группах
                    self.btnNorth.setChecked(True)
        if (9<=int(btnSenderName)<17): 
            self.lst_btn_group[1].setChecked(False)
            if (self.lst_btn_name[0+9].isChecked()&self.lst_btn_name[1+9].isChecked()&self.lst_btn_name[2+9].isChecked()&self.lst_btn_name[3+9].isChecked()&self.lst_btn_name[4+9].isChecked()&self.lst_btn_name[5+9].isChecked()&self.lst_btn_name[6+9].isChecked()&self.lst_btn_name[7+9].isChecked()):
                self.lst_btn_group[1].setChecked(True)
                if (self.lst_btn_group[0].isChecked()&self.lst_btn_group[1].isChecked()&self.lst_btn_group[2].isChecked()&self.lst_btn_group[3].isChecked()):# выделить север при всех нажатых группах
                    self.btnNorth.setChecked(True)
        if (17<=int(btnSenderName)<25): 
            self.lst_btn_group[2].setChecked(False)
            if (self.lst_btn_name[0+9+8].isChecked()&self.lst_btn_name[1+9+8].isChecked()&self.lst_btn_name[2+9+8].isChecked()&self.lst_btn_name[3+9+8].isChecked()&self.lst_btn_name[4+9+8].isChecked()&self.lst_btn_name[5+9+8].isChecked()&self.lst_btn_name[6+9+8].isChecked()&self.lst_btn_name[7+9+8].isChecked()):
                self.lst_btn_group[2].setChecked(True)
                if (self.lst_btn_group[0].isChecked()&self.lst_btn_group[1].isChecked()&self.lst_btn_group[2].isChecked()&self.lst_btn_group[3].isChecked()):# выделить север при всех нажатых группах
                    self.btnNorth.setChecked(True)
        if (25<=int(btnSenderName)<33): 
            self.lst_btn_group[3].setChecked(False)
            if (self.lst_btn_name[0+9+8+8].isChecked()&self.lst_btn_name[1+9+8+8].isChecked()&self.lst_btn_name[2+9+8+8].isChecked()&self.lst_btn_name[3+9+8+8].isChecked()&self.lst_btn_name[4+9+8+8].isChecked()&self.lst_btn_name[5+9+8+8].isChecked()&self.lst_btn_name[6+9+8+8].isChecked()&self.lst_btn_name[7+9+8+8].isChecked()):
                self.lst_btn_group[3].setChecked(True)
                if (self.lst_btn_group[0].isChecked()&self.lst_btn_group[1].isChecked()&self.lst_btn_group[2].isChecked()&self.lst_btn_group[3].isChecked()):# выделить север при всех нажатых группах
                    self.btnNorth.setChecked(True)
    def onClickbtnNorth(self): # Выбор антенн и групп
        if (self.btnNorth.isChecked()):
            for i in range(len(self.lst_btn_name)): self.lst_btn_name[i].setChecked(True)
            for i in range(len(self.lst_btn_group)): self.lst_btn_group[i].setChecked(True)
        else: 
            for i in range(len(self.lst_btn_name)): self.lst_btn_name[i].setChecked(False)
            for i in range(len(self.lst_btn_group)): self.lst_btn_group[i].setChecked(False)
    def onClickbtnGN01(self):
        if (self.lst_btn_group[0].isChecked()):
            for i in range(9): self.lst_btn_name[i].setChecked(True)
        else:
            for i in range(9): self.lst_btn_name[i].setChecked(False)
            self.btnNorth.setChecked(False)
        if (self.lst_btn_group[0].isChecked()&self.lst_btn_group[1].isChecked()&self.lst_btn_group[2].isChecked()&self.lst_btn_group[3].isChecked()):# выделить север при всех нажатых группах
            self.btnNorth.setChecked(True)    
    def onClickbtnGN02(self):
        if (self.lst_btn_group[1].isChecked()):
            for i in range(8): self.lst_btn_name[i+9].setChecked(True)
        else:
            for i in range(8): self.lst_btn_name[i+9].setChecked(False)
            self.btnNorth.setChecked(False)
        if (self.lst_btn_group[0].isChecked()&self.lst_btn_group[1].isChecked()&self.lst_btn_group[2].isChecked()&self.lst_btn_group[3].isChecked()):# выделить север при всех нажатых группах
            self.btnNorth.setChecked(True)
    def onClickbtnGN03(self):
        if (self.lst_btn_group[2].isChecked()):
            for i in range(8): self.lst_btn_name[i+9+8].setChecked(True)
        else:
            for i in range(8): self.lst_btn_name[i+9+8].setChecked(False)
            self.btnNorth.setChecked(False)
        if (self.lst_btn_group[0].isChecked()&self.lst_btn_group[1].isChecked()&self.lst_btn_group[2].isChecked()&self.lst_btn_group[3].isChecked()):# выделить север при всех нажатых группах
            self.btnNorth.setChecked(True)
    def onClickbtnGN04(self):
        if (self.lst_btn_group[3].isChecked()):
            for i in range(8): self.lst_btn_name[i+9+8+8].setChecked(True)
        else:
            for i in range(8): self.lst_btn_name[i+9+8+8].setChecked(False)
            self.btnNorth.setChecked(False)
        if (self.lst_btn_group[0].isChecked()&self.lst_btn_group[1].isChecked()&self.lst_btn_group[2].isChecked()&self.lst_btn_group[3].isChecked()):# выделить север при всех нажатых группах
            self.btnNorth.setChecked(True)
    def onClickbtnSunTracking(self):
        self.btnAzStop.setEnabled(False)# Disable buttons due to needs to stop SUN tracking
        self.btnAzLeft.setEnabled(False)
        self.btnAzRight.setEnabled(False)
        self.btnElStop.setEnabled(False)
        self.btnElUp.setEnabled(False)
        self.btnElDown.setEnabled(False)
        self.spinAzAngle.setDisabled(True);    
        for antName in dAntIP: 
            if (self.lst_btn_name[list(dAntIP.keys()).index(antName)].isChecked()): #если антенна выбрана
                self.setMode(dAntIP.get(antName),dSetAntMode.get('МСЭ-2'))# берем ее IP и посылыем туда команду
               
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
        for antName in dAntIP: 
            if (self.lst_btn_name[list(dAntIP.keys()).index(antName)].isChecked()): #если антенна выбрана
                self.setMode(dAntIP.get(antName),dSetAntMode.get('БСУ выключить'))# берем ее IP и посылыем туда команду
                
    def onClickbtnServicePoint(self):
        if (self.btnSunTracking.isChecked):
            self.btnSunTracking.setChecked(False)
        self.btnAzStop.setEnabled(False)
        self.btnAzLeft.setEnabled(False)
        self.btnAzRight.setEnabled(False)
        self.btnElStop.setEnabled(False)
        self.btnElUp.setEnabled(False)
        self.btnElDown.setEnabled(False)
        for antName in dAntIP: 
            if (self.lst_btn_name[list(dAntIP.keys()).index(antName)].isChecked()): #если антенна выбрана
                self.setAntCmd  (dAntIP.get(antName), OID_BSU + OID_setMode, Integer(dSetAntMode.get('БСУ выключить')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                self.setAntAngle(dAntIP.get(antName), OID_MM  + OID_setAngleAz, 0.1234) #Задаем угол по Азимуту
                self.setAntAngle(dAntIP.get(antName), OID_MM  + OID_setAngleEl, 0.6789) #Задаем угол по Углуместа
                self.setAntCmd  (dAntIP.get(antName), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                self.setAntCmd  (dAntIP.get(antName), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                #print(dAntIP.get(antName))
    def onClickAzForw(self):
        for antName in dAntIP: 
            if (self.lst_btn_name[list(dAntIP.keys()).index(antName)].isChecked()): #если антенна выбрана
                #setAntCmd  (dAntIP.get(antName), OID_BSU + OID_setMode, Integer(dSetAntMode.get('выключить')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                Angles = getAntAngle(dAntIP.get(antName))# Angles list: [0]- Azimuth, [1] - Elevation
                self.setAntAngle(dAntIP.get(antName), OID_MM  + OID_setAngleAz, Angles[0]+1.5) #Задаем угол по Азимуту
                self.setAntCmd  (dAntIP.get(antName), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду

    def onClickAzBack(self):
        for antName in dAntIP: 
            if (self.lst_btn_name[list(dAntIP.keys()).index(antName)].isChecked()): #если антенна выбрана
                #setAntCmd  (dAntIP.get(antName), OID_BSU + OID_setMode, Integer(dSetAntMode.get('выключить')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                Angles = getAntAngle(dAntIP.get(antName))# Angles list: [0]- Azimuth, [1] - Elevation
                self.setAntAngle(dAntIP.get(antName), OID_MM  + OID_setAngleAz, Angles[0]-1.5) #Задаем угол по Азимуту
                self.setAntCmd  (dAntIP.get(antName), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду

    def onClickAzStop(self):
        for antName in dAntIP: 
            if (self.lst_btn_name[list(dAntIP.keys()).index(antName)].isChecked()): #если антенна выбрана
                self.setAntCmd(dAntIP.get(antName), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Waiting')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду

    def onClickElUp(self):
        for antName in dAntIP: 
            if (self.lst_btn_name[list(dAntIP.keys()).index(antName)].isChecked()): #если антенна выбрана
                #setAntCmd  (dAntIP.get(antName), OID_BSU + OID_setMode, Integer(dSetAntMode.get('выключить')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                Angles = getAntAngle(dAntIP.get(antName))# Angles list: [0]- Azimuth, [1] - Elevation
                self.setAntAngle(dAntIP.get(antName), OID_MM  + OID_setAngleEl, Angles[1]+1.5) #Задаем угол по Азимуту
                self.setAntCmd  (dAntIP.get(antName), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду

    def onClickElDown(self):
        for antName in dAntIP: 
            if (self.lst_btn_name[list(dAntIP.keys()).index(antName)].isChecked()): #если антенна выбрана
                #setAntCmd  (dAntIP.get(antName), OID_BSU + OID_setMode, Integer(dSetAntMode.get('выключить')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                Angles = getAntAngle(dAntIP.get(antName))# Angles list: [0]- Azimuth, [1] - Elevation
                self.setAntAngle(dAntIP.get(antName), OID_MM  + OID_setAngleEl, Angles[1]-1.5) #Задаем угол по Азимуту
                self.setAntCmd  (dAntIP.get(antName), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду

    def onClickElStop(self):
        for antName in dAntIP: 
            if (self.lst_btn_name[list(dAntIP.keys()).index(antName)].isChecked()): #если антенна выбрана
                self.setAntCmd(dAntIP.get(antName), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Waiting')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
    
    def onClickAzSet(self):
        for antName in dAntIP: 
            if (self.lst_btn_name[list(dAntIP.keys()).index(antName)].isChecked()): #если антенна выбрана
                #setAntCmd  (dAntIP.get(antName), OID_BSU + OID_setMode, Integer(dSetAntMode.get('выключить')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                self.setAntAngle(dAntIP.get(antName), OID_MM  + OID_setAngleAz, self.spinAzAngle.value()+0.00000001) #Задаем угол по Азимуту
                self.setAntCmd  (dAntIP.get(antName), OID_MM  + OID_cmdSetAz, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
    
    def onClickElSet(self):
        for antName in dAntIP: 
            if (self.lst_btn_name[list(dAntIP.keys()).index(antName)].isChecked()): #если антенна выбрана
                #setAntCmd  (dAntIP.get(antName), OID_BSU + OID_setMode, Integer(dSetAntMode.get('выключить')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду
                self.setAntAngle(dAntIP.get(antName), OID_MM  + OID_setAngleEl, self.spinElAngle.value()+0.00000001) #Задаем угол по Азимуту
                self.setAntCmd  (dAntIP.get(antName), OID_MM  + OID_cmdSetEl, Integer(dSetCommand.get('Set Angle')))# берем ее IP и посылыем ТУДА(OID) соответствующюю команду

        
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    NAC = NorthAntCtrl()
    NAC.show()
    exit(app.exec_())
    #sys.exit(app.exec_())
    