# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 11:05:48 2016

@author: Alena
"""
#для корректного выхода из приложения.
import sys;
# для использования графического интерфейса, для работы с сетью.
from PyQt5 import QtGui,QtWidgets,QtCore,QtNetwork
import struct
from astropy.time import Time, TimeDelta
from datetime import datetime

ANT_FORWARD=0x17
ANT_BACKWARD=0x15
ANT_FAST_FORWARD=0x13
ANT_FAST_BACKWARD=0x11
ANT_STOP_X=0x31
ANT_UP=0x14
ANT_DOWN=0x16
ANT_FAST_UP=0x10
ANT_FAST_DOWN=0x12
ANT_STOP_Y=0x30

GROUP_FORWARD=0x27
GROUP_BACKWARD=0x25
GROUP_FAST_FORWARD=0x23
GROUP_FAST_BACKWARD=0x21
GROUP_STOP_X=0x41
GROUP_UP=0x24
GROUP_DOWN=0x26
GROUP_FAST_UP=0x20
GROUP_FAST_DOWN=0x22
GROUP_STOP_Y=0x40
BUK16_POWER_READ=104
BUK16_CTL_PORT=10124
BUK16_POWER_ON=120
BUK16_POWER_OFF=121
BUK16_SCHEDULE_ENABLE=115
BUK16_SCHEDULE_DISABLE=116
BUK16_SCHEDULE_UPDATE=114
BUK16_SCHEDULE_GET=119

class SrhAntennaControl(QtWidgets.QMainWindow):
    
    def closeEvent(self,event):
        close = QtWidgets.QMessageBox.question(self, 'SrhAntennaControl', "Are you sure?",QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
        if close == QtWidgets.QMessageBox.Yes:
            self.Timer.stop()

            timeStr=Time(datetime.utcnow(), scale='utc').iso.split('.')[0]
            a=timeStr.split('-')
            b=a[0]+a[1]+a[2]
            c=b.split(':')
            d='SrhAntennaControl_'+ c[0].replace(' ','_')+c[1]+c[2]+'.log'
            with open(d, 'w') as yourFile:
                yourFile.write(str(self.Log.toPlainText()))
                
            settings = QtCore.QSettings('SrhAntennaControl.conf',QtCore.QSettings.IniFormat)
            if (self.antPointing.isChecked()):
                settings.setValue('antPointing',True)
            else:
                settings.setValue('antPointing',False)
                
            if (self.toggleScheduleButton.isChecked()):
                settings.setValue('toggleSchedule',True)
            else:
                settings.setValue('toggleSchedule',False)
                
            if (self.hourModeTime.isChecked()):
                settings.setValue('hourMode',0)
            elif(self.hourModeAngle.isChecked()):
                settings.setValue('hourMode',1)
            elif(self.hourModeManual.isChecked()):
                settings.setValue('hourMode',2)
            elif(self.hourModeStopTime.isChecked()):
                settings.setValue('hourMode',3)
                
            if (self.declModeTime.isChecked()):
                settings.setValue('declMode',0)
            elif(self.declModeAngle.isChecked()):
                settings.setValue('declMode',1)
            elif(self.declModeManual.isChecked()):
                settings.setValue('declMode',2)
                
            settings.setValue('windowFrame',self.geometry())    
            event.accept()
        else:
            event.ignore()
        
        
    def onClickedButton(self,value):
        self.AntButton[self.antennaInArray].setChecked(False)
        antenna=int(self.sender().text())
        if antenna >= 49 and antenna <= 64:
            self.antennaInGroup=antenna-49
            self.antennaInArray=self.antennaInGroup
        elif  antenna >= 65 and antenna <= 80:
            self.antennaInGroup=antenna-65
            self.antennaInArray=self.antennaInGroup + 16
        elif  antenna >= 177 and antenna <= 192:
            self.antennaInGroup=antenna-177
            self.antennaInArray=self.antennaInGroup + 32
            
        self.AntButton[self.antennaInArray].setChecked(True)
        
           
    def onAntPointing(self, value):
        for i in range(48):
            self.AntButton[i].setVisible(value)
        if (value):            
            self.hourModeAngle.setChecked(True)
            self.declModeAngle.setChecked(True)
            
        self.hourModeTime.setEnabled(not value)         
        self.hourModeManual.setEnabled(not value)
        self.hourModeStopTime.setEnabled(not value)
        self.declModeTime.setEnabled(not value)
        self.declModeManual.setEnabled(not value)
        
        self.grButton4.setChecked(value)
        self.grButton5.setChecked(value)
        self.grButton12.setChecked(value)
        self.powerButton.setEnabled(not value)
        
        
        
    def onGrButton4(self, value):
        for i in range(16):
            self.AntButton[i].setChecked(0)
            
    def onGrButton5(self, value):
        for i in range(16):
            self.AntButton[i+16].setChecked(0)
            
    def onGrButton12(self, value):
        for i in range(16):
            self.AntButton[i+32].setChecked(0)
            
    def checkAntennaNumber(self,antennaNumber):
        result=False
        if(self.grButton4.isChecked() and self.AntNumber[self.antennaInArray] >= 49 and self.AntNumber[self.antennaInArray] <= 64):
            result=True
        elif(self.grButton5.isChecked() and self.AntNumber[self.antennaInArray] >= 65 and self.AntNumber[self.antennaInArray] <= 80):               
            result=True
        elif(self.grButton12.isChecked() and self.AntNumber[self.antennaInArray] >= 177 and self.AntNumber[self.antennaInArray] <= 192):               
            result=True
        return(result)
        
    def onFastForwForw(self):
        if (self.antPointing.isChecked()):
            if(self.checkAntennaNumber(self.antennaInArray)):
                self.hourInputAngleToTimer(7)
                self.pointingMode=7
                self.AntLabel[self.antennaInArray].setPixmap(self.antCmdIcons[ANT_FAST_FORWARD])
                self.sendToBUK(ANT_FAST_FORWARD)
        else:
            if (self.hourModeStopTime.isChecked()):            
                a=self.hourInput.toPlainText()
                b=a.split(":")
                stopSec=int(b[0])*3600.+int(b[1])*60.+int(b[2])
                curTime=Time.now().to_datetime()
                curSec=curTime.hour*3600+curTime.minute*60+curTime.second
                dt=(curSec-stopSec)/7
                self.MotionTimer.start(int(dt*1000))
                self.MotionTimeInterval=int(dt)
                self.pointingMode=7
                self.sendToBUK(GROUP_FAST_FORWARD)
            elif(self.hourModeAngle.isChecked()):
                self.hourInputAngleToTimer(7)
                self.pointingMode=7
                self.sendToBUK(GROUP_FAST_FORWARD)
            
    def onForw(self):
        if(self.hourModeTime.isChecked()):
           self.hourInputTimeToTimer()
        self.pointingMode=5
        self.sendToBUK(GROUP_FORWARD)        
        
    def onFastBack(self):        
        if(self.hourModeTime.isChecked()):
           self.hourInputTimeToTimer()
        self.pointingMode=2
        self.sendToBUK(GROUP_FAST_BACKWARD)
        
    def onFastForw(self):
        if(self.hourModeTime.isChecked()):
           self.hourInputTimeToTimer()
        self.pointingMode=6
        self.sendToBUK(GROUP_FAST_FORWARD)
        
    def onFastBackForw(self):
        if (self.antPointing.isChecked()):
            if(self.checkAntennaNumber(self.antennaInArray)):
                self.hourInputAngleToTimer(7)
                self.pointingMode=7
                self.AntLabel[self.antennaInArray].setPixmap(self.antCmdIcons[ANT_FAST_BACKWARD])
                self.sendToBUK(ANT_FAST_BACKWARD)
        else:
            if (self.hourModeStopTime.isChecked()):            
                a=self.hourInput.toPlainText()
                b=a.split(":")
                stopSec=int(b[0])*3600.+int(b[1])*60.+int(b[2])
                curTime=Time.now().to_datetime()
                curSec=curTime.hour*3600+curTime.minute*60+curTime.second
                dt=(curSec-stopSec)/9
                self.MotionTimer.start(int(dt*1000))
                self.MotionTimeInterval=int(dt)
                self.pointingMode=1
                self.sendToBUK(GROUP_FAST_BACKWARD)
            elif(self.hourModeAngle.isChecked()):
                self.hourInputAngleToTimer(9)
                self.pointingMode=1
                self.sendToBUK(GROUP_FAST_BACKWARD)
            
    def onBack(self):
        if(self.hourModeTime.isChecked()):
           self.hourInputTimeToTimer()
        self.pointingMode=3        
        self.sendToBUK(GROUP_BACKWARD)
    
    def declToTimer(self,speed):
       if(self.declModeTime.isChecked()):
           self.declInputTimeToTimer()
       elif(self.declModeAngle.isChecked()):
           self.declInputAngleToTimer(speed)        
            
    def onDown(self):
       self.pointingMode=10
       if (self.antPointing.isChecked()):
           if(self.checkAntennaNumber(self.antennaInArray)):
               self.declToTimer(1)           
               self.AntLabel[self.antennaInArray].setPixmap(self.antCmdIcons[ANT_DOWN])
               self.sendToBUK(ANT_DOWN)
       else:
           self.declToTimer()  
           self.sendToBUK(GROUP_DOWN)
        
    def onUp(self):
       self.pointingMode=9
       if (self.antPointing.isChecked()):
           if(self.checkAntennaNumber(self.antennaInArray)):
               self.declToTimer(1)           
               self.AntLabel[self.antennaInArray].setPixmap(self.antCmdIcons[ANT_UP])
               self.sendToBUK(ANT_UP)
       else:
           self.declToTimer()           
           self.sendToBUK(GROUP_UP)
        
    def onFastUp(self):        
       self.pointingMode=8
       if (self.antPointing.isChecked()):
           if(self.checkAntennaNumber(self.antennaInArray)):
               self.declToTimer(8)           
               self.AntLabel[self.antennaInArray].setPixmap(self.antCmdIcons[ANT_FAST_UP])
               self.sendToBUK(ANT_FAST_UP)
       else:
           self.declToTimer()           
           self.sendToBUK(GROUP_FAST_UP)
           
    def onFastDown(self):
       self.pointingMode=11
       if (self.antPointing.isChecked()):
           if(self.checkAntennaNumber(self.antennaInArray)):
               self.declToTimer(8)           
               self.AntLabel[self.antennaInArray].setPixmap(self.antCmdIcons[ANT_FAST_DOWN])
               self.sendToBUK(ANT_FAST_DOWN)
       else:
           self.declToTimer()           
           self.sendToBUK(GROUP_FAST_DOWN)
            
    def declInputTimeToTimer(self):
           a=self.declInput.toPlainText()
           b=a.split(":")
           sec=int(b[0])*3600.+int(b[1])*60.+int(b[2])
           self.MotionTimer.start(int(sec*1000))
           self.MotionTimeInterval=sec
           
    def declInputAngleToTimer(self,speed):
           a=self.declInput.toPlainText()
           b=a.split(":")
           sec=(int(b[0])*3600.+int(b[1])*60.+int(b[2]))/(15.*speed)
           self.MotionTimer.start(int(sec*1000))
           self.MotionTimeInterval=int(sec+0.5)
           
    def hourInputTimeToTimer(self):
           a=self.hourInput.toPlainText()
           b=a.split(":")
           sec=int(b[0])*3600.+int(b[1])*60.+int(b[2])
           self.MotionTimer.start(int(sec*1000))
           self.MotionTimeInterval=sec

           
    def hourInputAngleToTimer(self,speed):
           a=self.hourInput.toPlainText()
           b=a.split(":")
           sec=(int(b[0])*3600.+int(b[1])*60.+int(b[2]))/(15.*speed)
           self.MotionTimer.start(int(sec*1000))
           self.MotionTimeInterval=int(sec+0.5)
           
    def onStop(self):        
       self.pointingMode=4
       if (self.antPointing.isChecked()):
           self.AntLabel[self.antennaInArray].setPixmap(self.antCmdIcons[ANT_STOP_X])
           self.sendToBUK(ANT_STOP_X)
       else:
           self.sendToBUK(GROUP_STOP_X)
        
    def onStopDecl(self):        
       self.pointingMode=4
       if (self.antPointing.isChecked()):
           self.AntLabel[self.antennaInArray].setPixmap(self.antCmdIcons[ANT_STOP_Y])
           self.sendToBUK(ANT_STOP_Y)
       else:
           self.sendToBUK(GROUP_STOP_Y)
        
    def onHourModeTime(self,value):
        if(value):
            self.FastForwForw.setEnabled(False)
            self.FastBackForw.setEnabled(False)
            self.FastBack.setEnabled(True)
            self.FastForw.setEnabled(True)
            self.Forw.setEnabled(True)
            self.Back.setEnabled(True)
        
    def onHourModeStopTime(self,value):
        if(value):
            self.FastForwForw.setEnabled(True)
            self.FastBackForw.setEnabled(True)
            self.FastBack.setEnabled(False)
            self.FastForw.setEnabled(False)
            self.Forw.setEnabled(False)
            self.Back.setEnabled(False)
        
        
    def onHourModeAngle(self,value):
        if(value):
            self.FastForwForw.setEnabled(True)
            self.FastBackForw.setEnabled(True)
            self.FastBack.setEnabled(False)
            self.FastForw.setEnabled(False)
            self.Forw.setEnabled(False)
            self.Back.setEnabled(False)
        
    def onHourModeManual(self,value):
        if(value):
            self.FastForwForw.setEnabled(False)
            self.FastBackForw.setEnabled(False)
            self.FastBack.setEnabled(True)
            self.FastForw.setEnabled(True)
            self.Forw.setEnabled(True)
            self.Back.setEnabled(True)
        
    def onPowerButton(self):
        if(self.powerButton.isChecked()):
            self.powerButton.setStyleSheet('QPushButton {color: green}')
            self.sendToBUK(BUK16_POWER_ON)
        else:
            self.powerButton.setStyleSheet('QPushButton {color: gray}')
            self.sendToBUK(BUK16_POWER_OFF)
            
    def onUpdateScheduleButton(self):
        self.sendToBUK(BUK16_SCHEDULE_UPDATE)
        
    def onToggleScheduleButton(self,value):
        if(value):          
           self.sendToBUK(BUK16_SCHEDULE_ENABLE)
           self.toggleScheduleButton.setText('Enable')
        else:
           self.sendToBUK(BUK16_SCHEDULE_DISABLE)
           self.toggleScheduleButton.setText('Disable')
           
    def sendToBUK(self,cmd):
        timeStr=Time(datetime.utcnow(), scale='utc').iso.split('.')[0]
        msg = ('%s'%(self.BUKCommand[cmd])) 
        if (self.antPointing.isChecked()):
            self.antennaCurrentCmd=cmd
            antennaMask=bytearray(16)
            if (cmd==ANT_FAST_FORWARD or cmd==ANT_FAST_BACKWARD or cmd==ANT_BACKWARD or cmd==ANT_FORWARD or cmd==ANT_STOP_X):                
                antennaMask[self.antennaInGroup]=0x01            
            elif (cmd==ANT_DOWN or cmd==ANT_FAST_UP or cmd==ANT_FAST_DOWN or cmd==ANT_UP or cmd==ANT_STOP_Y):                
                antennaMask[self.antennaInGroup]=ANT_FAST_UP
            if (self.antennaInArray // 16 == 0):
                bukAddr=QtNetwork.QHostAddress('192.168.0.26')
                self.ClientBUK.writeDatagram(QtCore.QByteArray(struct.pack('<lll16s100s',cmd,128,0,antennaMask,bytes(100))),bukAddr,BUK16_CTL_PORT)
                self.Log.append(('%s, %s, Antenna_%d, %s'%(timeStr,self.IP_BUKToGroup[26],self.antennaInGroup + 49,msg)))
            elif (self.antennaInArray // 16 == 1):
                bukAddr=QtNetwork.QHostAddress('192.168.0.27')
                self.ClientBUK.writeDatagram(QtCore.QByteArray(struct.pack('<lll16s100s',cmd,128,0,antennaMask,bytes(100))),bukAddr,BUK16_CTL_PORT)
                self.Log.append(('%s, %s, Antenna_%d, %s'%(timeStr,self.IP_BUKToGroup[27],self.antennaInGroup + 65,msg)))
            elif (self.antennaInArray // 16 == 2):
                bukAddr=QtNetwork.QHostAddress('192.168.0.34')
                self.ClientBUK.writeDatagram(QtCore.QByteArray(struct.pack('<lll16s100s',cmd,128,0,antennaMask,bytes(100))),bukAddr,BUK16_CTL_PORT)
                self.Log.append(('%s, %s, Antenna_%d, %s'%(timeStr,self.IP_BUKToGroup[34],self.antennaInGroup + 177,msg)))
        else:
            if(self.grButton4.isChecked()):
                bukAddr=QtNetwork.QHostAddress('192.168.0.26')
                self.ClientBUK.writeDatagram(QtCore.QByteArray(struct.pack('<lllll108s',cmd,128,0,1,1,bytes(108))),bukAddr,BUK16_CTL_PORT)
                self.Log.append(('%s, %s, %s'%(timeStr,self.IP_BUKToGroup[26],msg)))
            if(self.grButton5.isChecked()):
                bukAddr=QtNetwork.QHostAddress('192.168.0.27')
                self.ClientBUK.writeDatagram(QtCore.QByteArray(struct.pack('<lllll108s',cmd,128,0,1,1,bytes(108))),bukAddr,BUK16_CTL_PORT)
                self.Log.append(('%s, %s, %s'%(timeStr,self.IP_BUKToGroup[27],msg)))
            if(self.grButton12.isChecked()):
                bukAddr=QtNetwork.QHostAddress('192.168.0.34')
                self.ClientBUK.writeDatagram(QtCore.QByteArray(struct.pack('<lllll108s',cmd,128,0,1,1,bytes(108))),bukAddr,BUK16_CTL_PORT)
                self.Log.append(('%s, %s, %s'%(timeStr,self.IP_BUKToGroup[34],msg)))

        
    def onMotionTimer(self):
        self.MotionTimer.stop()
        if (self.antPointing.isChecked()):
            if (self.antennaCurrentCmd==ANT_FAST_FORWARD or self.antennaCurrentCmd==ANT_FAST_BACKWARD or self.antennaCurrentCmd==ANT_BACKWARD):
                self.sendToBUK(ANT_FORWARD)                
            elif (self.antennaCurrentCmd==ANT_DOWN or self.antennaCurrentCmd==ANT_FAST_UP or self.antennaCurrentCmd==ANT_FAST_DOWN or self.antennaCurrentCmd==ANT_UP):
                self.sendToBUK(ANT_STOP_Y)
            self.AntLabel[self.antennaInArray].setPixmap(self.antCmdIcons[ANT_FORWARD])                
        else:
            if(self.pointingMode==1 or self.pointingMode==7):
                self.sendToBUK(GROUP_FORWARD)
            else:    
                self.sendToBUK(GROUP_STOP_Y)
        
    def onTimer(self):
        timeStr=Time(datetime.utcnow(), scale='utc').iso.split('.')[0]
        msg = ('%s'%(self.BUKCommand[BUK16_POWER_READ]))    
        bukAddr=QtNetwork.QHostAddress('192.168.0.26')
        self.ClientBUK.writeDatagram(QtCore.QByteArray(struct.pack('<lllll108s',BUK16_POWER_READ,128,0,1,1,bytes(108))),bukAddr,BUK16_CTL_PORT)
        self.Log.append(('%s, to %s, %s'%(timeStr,self.IP_BUKToGroup[26],msg)))
        bukAddr=QtNetwork.QHostAddress('192.168.0.27')
        self.ClientBUK.writeDatagram(QtCore.QByteArray(struct.pack('<lllll108s',BUK16_POWER_READ,128,0,1,1,bytes(108))),bukAddr,BUK16_CTL_PORT)
        self.Log.append(('%s, to %s, %s'%(timeStr,self.IP_BUKToGroup[27],msg)))
        bukAddr=QtNetwork.QHostAddress('192.168.0.34')
        self.ClientBUK.writeDatagram(QtCore.QByteArray(struct.pack('<lllll108s',BUK16_POWER_READ,128,0,1,1,bytes(108))),bukAddr,BUK16_CTL_PORT)        
        self.Log.append(('%s, to %s, %s'%(timeStr,self.IP_BUKToGroup[34],msg)))
        self.DateTime.setText(Time.now().utc.iso.split(".")[0])
        if(self.MotionTimeInterval > 0):
           print(self.MotionTimeInterval)
           self.MotionTimeInterval -= 1 
           deltaTime=TimeDelta(self.MotionTimeInterval, format='sec')
           remTime=Time('2001-01-01T00:00:00',format = 'isot')+deltaTime
           self.StopTime.setText(remTime.utc.iso.split(".")[0].split(" ")[1])
        else:
            curTime=Time.now()
            self.StopTime.setText(curTime.utc.iso.split(".")[0].split(" ")[1])
            
    def onScheduleTimer(self):
        bukAddr=QtNetwork.QHostAddress('192.168.0.26')
        self.ClientBUK.writeDatagram(QtCore.QByteArray(struct.pack('<lllll108s',BUK16_SCHEDULE_GET,128,0,1,1,bytes(108))),bukAddr,BUK16_CTL_PORT)

            
    def onClientBUKready(self):
        timeStr=Time(datetime.utcnow(), scale='utc').iso.split('.')[0]

        buf = self.ClientBUK.readDatagram(128)
        groupAddress = int(buf[1].toString().split('.')[3]);
        packet = struct.unpack('<l l l B B B B l l l l l l l l l l l l l l l l l l l l l l l L H H H H H H H H',QtCore.QByteArray(buf[0]));
        packetId=packet[0]
        groupCmd=packet[5]
        groupState=packet[3] & 0x15
        if(packetId==BUK16_POWER_READ):
            self.groupIcons[groupAddress].setPixmap(self.groupCmdIcons[groupCmd])
            self.groupStateIcons[groupAddress].setPixmap(self.groupStatePictures[groupState])
            if(groupState & 0x10):
                strPol = 'Polarization YES'
            else:
                strPol = 'Polarization NO'
                
            if(groupState & 0x05):
                strPow = 'Power NO'
            else:
                strPow = 'Power YES'
                
                
            msg = ('%s, %s, %s'%(self.BUKCommand[groupCmd],strPol,strPow))
            self.Log.append(('%s, from %s, %s'%(timeStr,self.IP_BUKToGroup[groupAddress],msg)))
        elif(packetId==BUK16_SCHEDULE_GET):
            packet = struct.unpack('<l l l l l l l l l l l l l l l l l l l l l l l l l l l l l l L H H',QtCore.QByteArray(buf[0]));
            if self.toggleScheduleButton.isChecked():
                self.scheduleText.setText(\
                Time(packet[3]+packet[4]*1e-6,format='unix',scale='utc').utc.iso.split()[1]+', '+\
                Time(packet[6]+packet[7]*1e-6,format='unix',scale='utc').utc.iso.split()[1]+', '+\
                Time(packet[9]+packet[10]*1e-6,format='unix',scale='utc').utc.iso.split()[1]+', '+\
                Time(packet[12]+packet[13]*1e-6,format='unix',scale='utc').utc.iso.split()[1]+', '+\
                Time(packet[15]+packet[16]*1e-6,format='unix',scale='utc').utc.iso.split()[1]+', '+\
                Time(packet[18]+packet[19]*1e-6,format='unix',scale='utc').utc.iso.split()[1]+', '+\
                Time(packet[21]+packet[22]*1e-6,format='unix',scale='utc').utc.iso.split()[1]+', '+\
                Time(packet[24]+packet[25]*1e-6,format='unix',scale='utc').utc.iso.split()[1])
            else:
                self.scheduleText.setText('Disabled')
                
        
    def __init__(self, parent=None):
        self.pointingMode=0
        self.antennaInGroup=0
        self.antennaInArray=0
        self.antennaCurrentCmd=0
        self.ephCurDate = Time(datetime.utcnow(), scale='utc')
        self.ephCurTime = Time(datetime.utcnow(), scale='utc')
        self.MotionTimeInterval=0
              
    # задает размер и положение окна.        
        QtWidgets.QMainWindow.__init__(self,parent);
#картинки       
        self.h_f=QtGui.QPixmap('ssrtGroups/res/h_f.png')#GROUP_FORWARD
        self.h_ff=QtGui.QPixmap('ssrtGroups/res/h_ff.png')#GROUP_FAST_FORWARD
        self.h_ff_f=QtGui.QPixmap('ssrtGroups/res/h_ff_f.png')
        self.h_fb_f=QtGui.QPixmap('ssrtGroups/res/h_fb_f.png')
        self.h_fb=QtGui.QPixmap('ssrtGroups/res/h_fb.png')#GROUP_FAST_BACKWARD
        self.h_f_f=QtGui.QPixmap('ssrtGroups/res/h_f_f.png')#GROUP_FAST_BACKWARD
        self.h_b=QtGui.QPixmap('ssrtGroups/res/h_b.png')#GROUP_BACKWARD
        self.hd_stop=QtGui.QPixmap('ssrtGroups/res/stop.png')#GROUP_STOP_X
        self.hd_stop1=QtGui.QPixmap('ssrtGroups/res/stop1.png')#GROUP_STOP_X
        self.d_u=QtGui.QPixmap('ssrtGroups/res/d_u.png')#GROUP_UP
        self.d_d=QtGui.QPixmap('ssrtGroups/res/d_d.png')#GROUP_DOWN
        self.d_fd=QtGui.QPixmap('ssrtGroups/res/d_fd.png')#GROUP_FAST_DOWN
        self.d_fu=QtGui.QPixmap('ssrtGroups/res/d_fu.png')#GROUP_FAST_UP
        self.s_nop=QtGui.QPixmap('ssrtGroups/res/no_p.png')#GROUP_FAST_UP
        self.s_no7p=QtGui.QPixmap('ssrtGroups/res/no_7p.png')#GROUP_FAST_UP
        self.s_no7=QtGui.QPixmap('ssrtGroups/res/no_7.png')#GROUP_FAST_UP
        self.s_ok=QtGui.QPixmap('ssrtGroups/res/ok.png')#GROUP_FAST_UP
 #словари   
        self.AntButton={}
        self.AntLabel={}
        self.AntNumber={}
        for i in range(16):            
           self.AntButton[i]=QtWidgets.QPushButton('%d'%(i+49),self)
           self.AntNumber[i]=i+49
           self.AntButton[i].setGeometry(110 + 40*i,50,30,16)
           self.AntButton[i].clicked.connect(self.onClickedButton)
           self.AntButton[i].setCheckable(True)
           self.AntButton[i].hide()
           self.AntLabel[i]=QtWidgets.QLabel(self)
           self.AntLabel[i].setPixmap(self.h_f)
           self.AntLabel[i].setGeometry(115 + 40*i,62,20,20)
           
           
           self.AntButton[i + 16]=QtWidgets.QPushButton('%d'%(i+65),self)
           self.AntNumber[i + 16]=i+65
           self.AntButton[i + 16].setGeometry(110 + 40*i,80,30,16)
           self.AntButton[i + 16].clicked.connect(self.onClickedButton)
           self.AntButton[i + 16].setCheckable(True)
           self.AntButton[i + 16].hide()
           self.AntLabel[i + 16]=QtWidgets.QLabel(self)
           self.AntLabel[i + 16].setPixmap(self.h_f)
           self.AntLabel[i + 16].setGeometry(115 + 40*i,92,20,20)
           
           self.AntButton[i + 32]=QtWidgets.QPushButton('%d'%(i+177),self)
           self.AntNumber[i + 32]=i+177
           self.AntButton[i + 32].setGeometry(110 + 40*i,110,30,16)
           self.AntButton[i + 32].clicked.connect(self.onClickedButton)
           self.AntButton[i + 32].setCheckable(True)
           self.AntButton[i + 32].hide()
           self.AntLabel[i + 32]=QtWidgets.QLabel(self)
           self.AntLabel[i + 32].setPixmap(self.h_f)
           self.AntLabel[i + 32].setGeometry(115 + 40*i,122,20,20)
           
        self.DateTime=QtWidgets.QLabel(self.ephCurDate.utc.iso.split(".")[0],self)
        self.DateTime.setGeometry(20,225,110,30)
        self.StopTime=QtWidgets.QLabel(self.ephCurTime.utc.iso.split(".")[0].split(" ")[1],self)
        self.StopTime.setGeometry(130,225,110,30)
        
        self.grLabel4=QtWidgets.QLabel(self)
        self.grLabel4.setGeometry(80,55,20,30)
        self.grState4=QtWidgets.QLabel(self)           
        self.grState4.setGeometry(80,52,20,20)
        self.grButton4=QtWidgets.QPushButton('Group4',self)
        self.grButton4.setGeometry(20,50,60,30)
        self.grButton4.setCheckable(True)
    
    
        self.grLabel5=QtWidgets.QLabel(self)
        self.grLabel5.setGeometry(80,85,20,30)
        self.grState5=QtWidgets.QLabel(self)
        self.grState5.setGeometry(80,82,20,20)
        self.grButton5=QtWidgets.QPushButton('Group5',self)
        self.grButton5.setGeometry(20,80,60,30)
        self.grButton5.setCheckable(True)
    
        self.grLabel12=QtWidgets.QLabel(self)
        self.grLabel12.setGeometry(80,115,20,30)
        self.grState12=QtWidgets.QLabel(self)
        self.grState12.setGeometry(80,112,20,20)
        self.grButton12=QtWidgets.QPushButton('Group12',self)
        self.grButton12.setGeometry(20,110,60,30)
        self.grButton12.setCheckable(True)
        
        self.powerLabel=QtWidgets.QLabel(self)
        self.powerLabel.setGeometry(80,115,20,30)
        self.powerState=QtWidgets.QLabel(self)
        self.powerState.setGeometry(80,112,20,20)
        self.powerButton=QtWidgets.QPushButton('ON/OFF',self)
        self.powerButton.setGeometry(730,10,60,25)
        self.powerButton.setCheckable(True)        
        self.powerButton.clicked.connect(self.onPowerButton)
        pal = self.powerButton.palette()
        pal.setColor(QtGui.QPalette.Button,QtGui.QColor('red'))
        self.powerButton.setAutoFillBackground(True)
        self.powerButton.setStyleSheet('QPushButton {color: green}') 
        self.powerButton.setPalette(pal)
        self.powerButton.update()
        
        self.toggleScheduleButton=QtWidgets.QPushButton('Disable',self)
        self.toggleScheduleButton.setGeometry(680,10,50,25)
        self.toggleScheduleButton.setCheckable(True)
        self.toggleScheduleButton.toggled.connect(self.onToggleScheduleButton)
        self.updateScheduleButton=QtWidgets.QPushButton('Update',self)        
        self.updateScheduleButton.setGeometry(630,10,50,25)
        self.updateScheduleButton.clicked.connect(self.onUpdateScheduleButton)
        
        self.hourModeTime=QtWidgets.QRadioButton('time',self)
        self.hourModeTime.toggled.connect(self.onHourModeTime)
        self.hourModeAngle=QtWidgets.QRadioButton('angle',self)
        self.hourModeAngle.toggled.connect(self.onHourModeAngle)
        self.hourModeManual=QtWidgets.QRadioButton('manual',self)
        self.hourModeManual.toggled.connect(self.onHourModeManual)
        self.hourModeStopTime=QtWidgets.QRadioButton('stop time',self)
        self.hourModeStopTime.toggled.connect(self.onHourModeStopTime)
        self.hourModeGroup=QtWidgets.QButtonGroup(self)
        self.hourModeGroup.addButton(self.hourModeTime) 
        self.hourModeGroup.addButton(self.hourModeAngle) 
        self.hourModeGroup.addButton(self.hourModeManual) 
        self.hourModeGroup.addButton(self.hourModeStopTime)
        self.hourModeTime.setGeometry(110,160,50,30) 
        self.hourModeAngle.setGeometry(170,160,50,30) 
        self.hourModeManual.setGeometry(230,160,60,30) 
        self.hourModeStopTime.setGeometry(290,160,70,30)
        self.hourInput=QtWidgets.QTextEdit(self);        
        self.hourInput.setGeometry(110,190,100,30); 
        self.hourInput.setText('00:30:00')
        self.scheduleText=QtWidgets.QTextEdit(self);        
        self.scheduleText.setGeometry(10,10,620,25); 
       
       
        self.declModeTime=QtWidgets.QRadioButton('time',self)
        self.declModeAngle=QtWidgets.QRadioButton('angle',self)
        self.declModeManual=QtWidgets.QRadioButton('manual',self)
        self.declModeGroup=QtWidgets.QButtonGroup(self)
        self.declModeGroup.addButton(self.declModeTime) 
        self.declModeGroup.addButton(self.declModeAngle) 
        self.declModeGroup.addButton(self.declModeManual) 
        self.declModeTime.setGeometry(400,160,50,30) 
        self.declModeAngle.setGeometry(470,160,50,30) 
        self.declModeManual.setGeometry(540,160,60,30)
        self.declInput=QtWidgets.QTextEdit(self); 
        self.declInput.setGeometry(400,190,100,30); 
        self.declInput.setText('00:30:00')
       
        self.groupStateIcons={}
        self.groupStateIcons[22]=self.grState4
        self.groupStateIcons[26]=self.grState4        
        self.groupStateIcons[27]=self.grState5
        self.groupStateIcons[34]=self.grState12
        
        self.IP_BUKToGroup={}
        self.IP_BUKToGroup[22]='TestGroup'
        self.IP_BUKToGroup[26]='Group_04'
        self.IP_BUKToGroup[27]='Group_05'
        self.IP_BUKToGroup[34]='Group_12'
        
        self.BUKCommand={}
        self.BUKCommand[ANT_FORWARD]='ANT_FORWARD'
        self.BUKCommand[ANT_BACKWARD]='ANT_BACKWARD'
        self.BUKCommand[ANT_FAST_FORWARD]='ANT_FAST_FORWARD'      
        self.BUKCommand[ANT_FAST_BACKWARD]='ANT_FAST_BACKWARD'
        self.BUKCommand[ANT_STOP_X]='ANT_STOP_X'
        self.BUKCommand[ANT_UP]='ANT_UP'
        self.BUKCommand[ANT_DOWN]='ANT_DOWN'
        self.BUKCommand[ANT_FAST_UP]='ANT_FAST_UP'
        self.BUKCommand[ANT_FAST_DOWN]='ANT_FAST_DOWN'
        self.BUKCommand[ANT_STOP_Y]='ANT_STOP_Y'
        
        self.BUKCommand[GROUP_FORWARD]='GROUP_FORWARD'
        self.BUKCommand[GROUP_BACKWARD]='GROUP_BACKWARD'
        self.BUKCommand[GROUP_FAST_FORWARD]='GROUP_FAST_FORWARD'      
        self.BUKCommand[GROUP_FAST_BACKWARD]='GROUP_FAST_BACKWARD'
        self.BUKCommand[GROUP_STOP_X]='GROUP_STOP_X'
        self.BUKCommand[GROUP_UP]='GROUP_UP'
        self.BUKCommand[GROUP_DOWN]='GROUP_DOWN'
        self.BUKCommand[GROUP_FAST_UP]='GROUP_FAST_UP'
        self.BUKCommand[GROUP_FAST_DOWN]='GROUP_FAST_DOWN'
        self.BUKCommand[GROUP_STOP_Y]='GROUP_STOP_Y'
        
        self.BUKCommand[BUK16_POWER_READ]='BUK16_POWER_READ'
        self.BUKCommand[BUK16_POWER_ON]='BUK16_POWER_ON'
        self.BUKCommand[BUK16_POWER_OFF]='BUK16_POWER_OFF'
        self.BUKCommand[BUK16_SCHEDULE_ENABLE]='BUK16_SCHEDULE_ENABLE'
        self.BUKCommand[BUK16_SCHEDULE_DISABLE]='BUK16_SCHEDULE_DISABLE'
        self.BUKCommand[BUK16_SCHEDULE_UPDATE]='BUK16_SCHEDULE_UPDATE'
        
        

        self.groupStatePictures={}
        self.groupStatePictures[ANT_BACKWARD]=self.s_nop
        self.groupStatePictures[0x05]=self.s_no7p
        self.groupStatePictures[0x00]=self.s_no7
        self.groupStatePictures[ANT_FAST_UP]=self.s_ok
       
        self.groupIcons={}
        self.groupIcons[22]=self.grLabel4
        self.groupIcons[26]=self.grLabel4
        self.groupIcons[27]=self.grLabel5
        self.groupIcons[34]=self.grLabel12
       
        self.groupCmdIcons={}
        self.groupCmdIcons[GROUP_FORWARD]=self.h_f
        self.groupCmdIcons[GROUP_BACKWARD]=self.h_b
        self.groupCmdIcons[GROUP_FAST_FORWARD]=self.h_ff
        self.groupCmdIcons[GROUP_FAST_BACKWARD]=self.h_fb           
        self.groupCmdIcons[GROUP_UP]=self.d_u
        self.groupCmdIcons[GROUP_DOWN]=self.d_d
        self.groupCmdIcons[GROUP_FAST_UP]=self.d_fu
        self.groupCmdIcons[GROUP_FAST_DOWN]=self.d_fd
        self.groupCmdIcons[GROUP_STOP_Y]=self.hd_stop
        self.groupCmdIcons[GROUP_STOP_X]=self.hd_stop
        
        self.antCmdIcons={}
        self.antCmdIcons[ANT_FORWARD]=self.h_f
        self.antCmdIcons[ANT_BACKWARD]=self.h_b
        self.antCmdIcons[ANT_FAST_FORWARD]=self.h_ff
        self.antCmdIcons[ANT_FAST_BACKWARD]=self.h_fb           
        self.antCmdIcons[ANT_UP]=self.d_u
        self.antCmdIcons[ANT_DOWN]=self.d_d
        self.antCmdIcons[ANT_FAST_UP]=self.d_fu
        self.antCmdIcons[ANT_FAST_DOWN]=self.d_fd
        self.antCmdIcons[ANT_STOP_Y]=self.hd_stop
        self.antCmdIcons[ANT_STOP_X]=self.hd_stop
                    
        self.Timer=QtCore.QTimer(self)
        self.Timer.timeout.connect(self.onTimer)
        self.MotionTimer=QtCore.QTimer(self)
        self.MotionTimer.timeout.connect(self.onMotionTimer)
        
        self.ScheduleTimer=QtCore.QTimer(self)
        self.ScheduleTimer.timeout.connect(self.onScheduleTimer)
        
        self.ClientBUK=QtNetwork.QUdpSocket()
        self.ClientBUK.bind(QtNetwork.QHostAddress.Any,0)
        self.ClientBUK.readyRead.connect(self.onClientBUKready)
        
        self.antPointing=QtWidgets.QCheckBox("Pointing",self)
        self.antPointing.setGeometry(20,160,60,30)        
        self.antPointing.stateChanged.connect(self.onAntPointing)
        
        self.FastForwForw=QtWidgets.QPushButton(self)
        self.FastForwForw.setIcon(QtGui.QIcon(self.h_ff_f))
        self.FastForwForw.setIconSize(QtCore.QSize(20,20))
        self.FastForwForw.setGeometry(290,140,25,20)
        self.FastForwForw.clicked.connect(self.onFastForwForw)
        
        self.Forw=QtWidgets.QPushButton(self)
        self.Forw.setIcon(QtGui.QIcon(self.h_f_f))
        self.Forw.setIconSize(QtCore.QSize(20,50))
        self.Forw.setGeometry(230,140,25,20)
        self.Forw.clicked.connect(self.onForw)        
        
        self.FastBackForw=QtWidgets.QPushButton(self)
        self.FastBackForw.setIcon(QtGui.QIcon(self.h_fb_f))
        self.FastBackForw.setIconSize(QtCore.QSize(20,20))
        self.FastBackForw.setGeometry(110,140,25,20)
        self.FastBackForw.clicked.connect(self.onFastBackForw)
        
        self.Back=QtWidgets.QPushButton(self)
        self.Back.setIcon(QtGui.QIcon(self.h_b))
        self.Back.setIconSize(QtCore.QSize(20,20))
        self.Back.setGeometry(170,140,25,20)
        self.Back.clicked.connect(self.onBack)
        
        self.FastBack=QtWidgets.QPushButton(self)
        self.FastBack.setIcon(QtGui.QIcon(self.h_fb))
        self.FastBack.setIconSize(QtCore.QSize(20,20))
        self.FastBack.setGeometry(140,140,25,20)
        self.FastBack.clicked.connect(self.onFastBack)
        
        self.FastForw=QtWidgets.QPushButton(self)
        self.FastForw.setIcon(QtGui.QIcon(self.h_ff))
        self.FastForw.setIconSize(QtCore.QSize(20,20))
        self.FastForw.setGeometry(260,140,25,20)
        self.FastForw.clicked.connect(self.onFastForw)
        
        self.FastDown=QtWidgets.QPushButton(self)
        self.FastDown.setIcon(QtGui.QIcon(self.d_fd))
        self.FastDown.setIconSize(QtCore.QSize(20,50))
        self.FastDown.setGeometry(520,140,25,20)
        self.FastDown.clicked.connect(self.onFastDown)
        
        self.Down=QtWidgets.QPushButton(self)
        self.Down.setIcon(QtGui.QIcon(self.d_d))
        self.Down.setIconSize(QtCore.QSize(20,50))
        self.Down.setGeometry(490,140,25,20)
        self.Down.clicked.connect(self.onDown)
        
        self.FastUp=QtWidgets.QPushButton(self)
        self.FastUp.setIcon(QtGui.QIcon(self.d_fu))
        self.FastUp.setIconSize(QtCore.QSize(20,20))
        self.FastUp.setGeometry(400,140,25,20)
        self.FastUp.clicked.connect(self.onFastUp)
        
        self.Up=QtWidgets.QPushButton(self)
        self.Up.setIcon(QtGui.QIcon(self.d_u))
        self.Up.setIconSize(QtCore.QSize(20,20))
        self.Up.setGeometry(430,140,25,20)
        self.Up.clicked.connect(self.onUp)
        
        self.Stop=QtWidgets.QPushButton(self)
        self.Stop.setIcon(QtGui.QIcon(self.hd_stop1))
        self.Stop.setIconSize(QtCore.QSize(20,20))
        self.Stop.setGeometry(200,140,25,20)
        self.Stop.clicked.connect(self.onStop)
        
        self.StopDecl=QtWidgets.QPushButton(self)
        self.StopDecl.setIcon(QtGui.QIcon(self.hd_stop1))
        self.StopDecl.setIconSize(QtCore.QSize(20,20))
        self.StopDecl.setGeometry(460,140,25,20)
        self.StopDecl.clicked.connect(self.onStopDecl)
        

#подключение сигнала        
        self.grButton4.toggled.connect(self.onGrButton4)
        self.grButton5.toggled.connect(self.onGrButton5)
        self.grButton12.toggled.connect(self.onGrButton12)
        
#окно для отображения ввенных и полученных сообщений;    
        self.Log=QtWidgets.QTextEdit(self); 
        self.Log.setGeometry(20,250,700,90);
        
        settings = QtCore.QSettings('SrhAntennaControl.conf',QtCore.QSettings.IniFormat)
        
        if(settings.value('antPointing')=='true'):
            self.antPointing.setChecked(True)
            
        if(settings.value('toggleSchedule')=='true'):
            self.toggleScheduleButton.setChecked(True)            
            
        if(settings.value('hourMode')=='0'):
            self.hourModeTime.setChecked(True)
        elif(settings.value('hourMode')=='1'):            
            self.hourModeAngle.setChecked(True)
        elif(settings.value('hourMode')=='2'):
            self.hourModeManual.setChecked(True)
        elif(settings.value('hourMode')=='3'):
            self.hourModeStopTime.setChecked(True)
            
        if(settings.value('declMode')=='0'):
            self.declModeTime.setChecked(True)
        elif(settings.value('declMode')=='1'):
            self.declModeAngle.setChecked(True)
        elif(settings.value('declMode')=='2'):
            self.declModeManual.setChecked(True)
            
        self.setGeometry(settings.value('windowFrame'))

        self.Timer.start(900)
        self.ScheduleTimer.start(1000)
#завершение описания класса            
# для создания приложений и его многократного запуска.  
application = QtWidgets.QApplication.instance();
if not application:
    application = QtWidgets.QApplication(sys.argv);
    
if sys.platform == 'linux':
    font = QtGui.QFont()
    application.setFont(QtGui.QFont(font.defaultFamily(),8));
SRHantCtrl=SrhAntennaControl();
SRHantCtrl.setWindowTitle('SRH antenna control')
SRHantCtrl.show();
sys.exit(application.exec_());


