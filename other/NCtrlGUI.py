# -*- coding: utf-8 -*-
"""
Created on Fri May 29 23:25:15 2020

@author: User
"""



import sys
from PyQt5 import QtWidgets#QMainWindow, QAction, qApp, QApplication
from PyQt5.QtGui import QIcon

StyleSheet = '''
QPushButton::hover{background-color : #E1E1E1;}
QPushButton::checked{background-color :  solid rgb(100,150,180)}
'''

            
def initUI(self):
    exitAct = QtWidgets.QAction(QIcon('exit.png'), '&Exit', self)
    exitAct.setShortcut('Ctrl+Q')
    exitAct.setStatusTip('Exit application')
    exitAct.triggered.connect(self.close)
    #exitAct.triggered.connect(QtWidgets.qApp.quit)
    
    self.setWindowTitle('SRH 3-6 antenna control')
    self.statusBar()
    menubar = self.menuBar()
    menuFile = menubar.addMenu('&File')
    menuPref = menubar.addMenu('&Preferences')
    menuAbout = menubar.addMenu('&About')
    menuFile.addAction(exitAct)
    
    self.mainWidget = QtWidgets.QWidget() #Заводим основной виджет
    self.setCentralWidget(self.mainWidget) # и паркуем его на главном окне
    
    #self.list_of_names      = ['C00'] +['N0'+str(i) for i in range(1,33)] # Создание списка номеров антенн
    self.list_of_ant_North  = ['C00'] +['N0'+str(i) for i in range(1,33)] # Создание списка номеров антенн
    self.list_of_ant_West   = ['W0' + str(i) for i in range(1,10)]+['W' + str(i) for i in range(10,33)]# Создание списка номеров антенн по Западу
    self.list_of_ant_East   = ['E0' + str(i) for i in range(1,10)]+['E' + str(i) for i in range(10,33)] # Создание списка номеров антенн по Ближнему Востоку
    self.list_of_ant_FarEast= ['E'  + str(i) for i in range(33,65)]# Создание списка номеров антенн по Дальнему Востоку
    self.list_of_directions = ['North','West', 'East', 'FarEast'] # Создание списка имен плечей СРГ
 
    self.btnCentral = QtWidgets.QPushButton("C00")
    #self.btn3 = QtWidgets.QPushButton("Aramis")
    self.btnSRH = QtWidgets.QPushButton("SRH 3-6")
    self.btnSRH.setMinimumWidth(35*33 )#Устанавливаем минимальную ширину кнопок
    self.btnSRH.setCheckable(True)
    self.btnSRH.clicked.connect(self.onClickbtnSRH)
    self.btnSRH.setStyleSheet(StyleSheet)
 
    self.vbox         = QtWidgets.QVBoxLayout(self.mainWidget, spacing = 0)# Главный вертикальный бокс
    self.hbox_antenna = QtWidgets.QHBoxLayout(spacing = 0)
    self.vbox_antenna = QtWidgets.QVBoxLayout(spacing = 0)
    self.hbox_direction=QtWidgets.QHBoxLayout(spacing = 0)
    self.hbox_control = QtWidgets.QHBoxLayout(spacing = 0)
    self.grid_azimuth = QtWidgets.QGridLayout(spacing = 0) #Создаем widget сетку
    self.grid_elevation=QtWidgets.QGridLayout(spacing = 0) #Создаем widget сетку
    self.grid_mode    = QtWidgets.QGridLayout(spacing = 0) #Создаем widget сетку
    self.grid_log     = QtWidgets.QVBoxLayout(spacing = 0) #Создаем widget сетку
    #self.hbox_antenna.addStretch(1) #убираем промежуток между кнопками
    #Frame = QtWidgets.QFrame(w)
    self.antennaBox   = QtWidgets.QGroupBox('Antennas')#Группируем антенны в рамку 
    self.controlBox   = QtWidgets.QGroupBox('Control') #Группируем в рамку элементы управления 
    self.modeBox      = QtWidgets.QGroupBox('Mode')#Создаем рамку для  Режима
    self.azimuthBox   = QtWidgets.QGroupBox('Azimuth')#Создаем рамку для  Азимута
    self.elevationBox = QtWidgets.QGroupBox('Elevation')#Создаем рамку для  Угла места
    self.logBox       = QtWidgets.QGroupBox('log')#Создаем рамку для  Состояния
    self.vbox.addWidget(self.antennaBox)
    self.vbox.addLayout(self.hbox_control)
    
    self.hbox_control.addWidget(self.modeBox)
    self.hbox_control.addWidget(self.azimuthBox)
    self.hbox_control.addWidget(self.elevationBox)
    self.hbox_control.addWidget(self.logBox)
    
    self.azimuthBox.setLayout(self.grid_azimuth)# размещаем сетку в рамке Азимута
    self.elevationBox.setLayout(self.grid_elevation)# размещаем сетку в рамке Азимута
    self.modeBox.setLayout(self.grid_mode)# размещаем сетку в рамке Режима
    self.logBox.setLayout(self.grid_log)# размещаем сетку в рамке Режима
       
    self.btnSunTracking = QtWidgets.QPushButton("Sun tracking")
#    self.btnSunTracking.setCheckable(True)
    self.btnSunTracking.setStyleSheet(StyleSheet)#Устанавливаем фирменный стиль)
    self.grid_mode.addWidget(self.btnSunTracking,1,1)# Добавляем кнопку в сетку с коорд 1-1
    self.btnSunTracking.clicked.connect(self.onClickbtnSunTracking)# привязываем сигнал нажатия к ф-ии onClickbtnSunTracking
    
    self.btnStopMoving = QtWidgets.QPushButton("Stop moving")
    self.btnStopMoving.setCheckable(False)
    self.btnStopMoving.setStyleSheet(StyleSheet)#Устанавливаем фирменный стиль)
    self.grid_mode.addWidget(self.btnStopMoving,2,1)# Добавляем кнопку в сетку с коорд 1-1
    self.btnStopMoving.clicked.connect(self.onClickbtnStopMoving)# привязываем сигнал нажатия к ф-ии onClickbtnSunTracking
    
    self.btnServicePoint = QtWidgets.QPushButton("Service point")
    self.btnServicePoint.setCheckable(False)
    self.btnServicePoint.setStyleSheet(StyleSheet)#Устанавливаем фирменный стиль)
    self.grid_mode.addWidget(self.btnServicePoint,3,1)# Добавляем кнопку в сетку с коорд 1-1
    self.btnServicePoint.clicked.connect(self.onClickbtnServicePoint)# привязываем сигнал нажатия к ф-ии onClickbtnSunTracking
    
    self.btnAzRight = QtWidgets.QPushButton("Forward 1.5\N{DEGREE SIGN}")
    self.btnAzRight.setCheckable(False)
    self.btnAzRight.clicked.connect(self.onClickAzForw)
    self.btnAzRight.setStyleSheet(StyleSheet)#Устанавливаем фирменный стиль)
    self.grid_azimuth.addWidget(self.btnAzRight,3,3)# Добавляем кнопку в сетку с коорд 1-1

    self.btnAzStop = QtWidgets.QPushButton("Stop")
    self.btnAzStop.setCheckable(False)
    self.btnAzStop.clicked.connect(self.onClickAzStop)
    self.btnAzStop.setStyleSheet(StyleSheet)
    self.grid_azimuth.addWidget(self.btnAzStop,3,2)
    
    self.btnAzLeft = QtWidgets.QPushButton("Backward 1.5\N{DEGREE SIGN}")
    self.btnAzLeft.setCheckable(False)
    self.btnAzLeft.clicked.connect(self.onClickAzBack)
    self.btnAzLeft.setStyleSheet(StyleSheet)
    self.grid_azimuth.addWidget(self.btnAzLeft,3,1)
    
    self.spinAzAngle = QtWidgets.QDoubleSpinBox()
    self.spinAzAngle.setWrapping(True)
    #sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred)
    #self.lineAzAngle.setSizePolicy(sizePolicy)
    #self.spinAzAngle.setFixedWidth(80)
    self.grid_azimuth.addWidget(self.spinAzAngle,1,1)# Добавляем спинбокс в сетку с коорд 1-1
    self.spinAzAngle.setRange(-180.5,180.5)
    self.spinAzAngle.setSingleStep(0.05)
    self.spinAzAngle.setDecimals(3)#Количество знаков после запятой
    self.spinAzAngle.lineEdit().deselect()
    
    self.btnAzSet = QtWidgets.QPushButton("Set Azimuth")
    self.btnAzSet.clicked.connect(self.onClickAzSet)
    self.btnAzSet.setStyleSheet(StyleSheet)
    self.grid_azimuth.addWidget(self.btnAzSet,1,2)
    
    self.spinElAngle = QtWidgets.QDoubleSpinBox()
    self.spinElAngle.setWrapping(True)
    self.spinElAngle.setFixedWidth(80)
    self.grid_elevation.addWidget(self.spinElAngle,1,1)# Добавляем спинбокс в сетку с коорд 1-1
    self.spinElAngle.setRange(-0.5,90.5)
    self.spinElAngle.setSingleStep(0.05)
    self.spinElAngle.setDecimals(3)#Количество знаков после запятой
    self.spinElAngle.lineEdit().deselect()
    
    self.btnElSet = QtWidgets.QPushButton("Set Elevation")
    self.btnElSet.clicked.connect(self.onClickElSet)
    self.btnElSet.setStyleSheet(StyleSheet)
    self.grid_elevation.addWidget(self.btnElSet,2,1)
    
    self.btnElUp = QtWidgets.QPushButton("Up 1.5\N{DEGREE SIGN}")
    self.btnElUp.setCheckable(False)
    self.btnElUp.clicked.connect(self.onClickElUp)
    self.btnElUp.setStyleSheet(StyleSheet)#Устанавливаем фирменный стиль)
    self.btnElUp.setFixedWidth(80)
    self.grid_elevation.addWidget(self.btnElUp,1,3)# Добавляем кнопку в сетку с коорд 1-1

    self.btnElStop = QtWidgets.QPushButton("Stop")
    self.btnElStop.setCheckable(False)
    self.btnElStop.clicked.connect(self.onClickElStop)
    self.btnElStop.setStyleSheet(StyleSheet)
    self.btnElStop.setFixedWidth(80)
    self.grid_elevation.addWidget(self.btnElStop,2,3)
    
    self.btnElDown = QtWidgets.QPushButton("Down 1.5\N{DEGREE SIGN}")
    self.btnElDown.setCheckable(False)
    self.btnElDown.clicked.connect(self.onClickElDown)
    self.btnElDown.setStyleSheet(StyleSheet)
    self.btnElDown.setFixedWidth(80)
    self.grid_elevation.addWidget(self.btnElDown,3,3)
    
    self.log = QtWidgets.QTextEdit()
    self.grid_log.addWidget(self.log)
    sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred)
    self.log.setFixedHeight(400)
    self.log.setSizePolicy(sizePolicy)
    self.log.setFixedHeight(80)
    self.log.setReadOnly(True)
    #self.grid_control.addWidget(self.btn3,0,2,0,-1)
    
    #self.grid_control.addWidget(self.btn1,5,0)
    #print(self.grid_control.columnCount())
    self.antennaBox.setLayout(self.vbox_antenna)   #Вкладываем вертикальный бокс
    self.vbox_antenna.addWidget(self.btnSRH)     # в него кнопку выбора луча
    self.vbox_antenna.addLayout(self.hbox_direction)   # затем горизонтальный для плечей
    self.vbox_antenna.addLayout(self.hbox_antenna) # затем горизонтальный для антенн
    self.lst_btn_name = [] # Создание списка кнопок номеров антенн
    self.lst_btn_ant_North  = [] # Создание списка кнопок номеров антенн Севера
    self.lst_btn_ant_West   = [] # Создание списка кнопок номеров антенн Запада
    self.lst_btn_ant_East   = [] # Создание списка кнопок номеров антенн Востока
    self.lst_btn_ant_FarEast= [] # Создание списка кнопок номеров антенн Дальнего Востока
    
#    for name in self.list_of_names: # Генерация кнопок по списку антенн
#            btn = QtWidgets.QPushButton(name)
#            btn.setMinimumWidth(35 )#Устанавливаем минимальную ширину кнопок
#            btn.setCheckable(True)
#            btn.setVisible(False)
#            #btn.clicked.connect(self.pr)
#            self.hbox_antenna.addWidget(btn)
#            btn.setStyleSheet(StyleSheet)
#            btn.setObjectName(str(self.list_of_names.index(name)))# устанавливавем имя обьекта равное его индексу в списке
#            btn.clicked.connect(self.onClickbtnName)
#            self.lst_btn_name.append(btn) # Формируем список кнопок номеров антенн
    
    for name in self.list_of_ant_North: # Генерация кнопок по списку антенн
            btn = QtWidgets.QPushButton(name)
            btn.setMinimumWidth(35 )#Устанавливаем минимальную ширину кнопок
            btn.setCheckable(True)
            btn.setVisible(False)
            #btn.clicked.connect(self.pr)
            self.hbox_antenna.addWidget(btn)
            btn.setStyleSheet(StyleSheet)
            btn.setObjectName(str(self.list_of_ant_North.index(name)))# устанавливавем имя обьекта равное его индексу в списке
            #btn.clicked.connect(self.onClickbtnAntNorth)
            self.lst_btn_ant_North.append(btn) # Формируем список кнопок номеров антенн
    
    for name in self.list_of_ant_West: # Генерация кнопок по списку антенн
            btn = QtWidgets.QPushButton(name)
            btn.setMinimumWidth(35 )#Устанавливаем минимальную ширину кнопок
            btn.setCheckable(True)
            btn.setVisible(False)
            #btn.clicked.connect(self.pr)
            self.hbox_antenna.addWidget(btn)
            btn.setStyleSheet(StyleSheet)
            btn.setObjectName(str(self.list_of_ant_West.index(name)))# устанавливавем имя обьекта равное его индексу в списке
            #btn.clicked.connect(self.onClickbtnAntWest)
            self.lst_btn_ant_West.append(btn) # Формируем список кнопок номеров антенн
    for name in self.list_of_ant_East: # Генерация кнопок по списку антенн
            btn = QtWidgets.QPushButton(name)
            btn.setMinimumWidth(35 )#Устанавливаем минимальную ширину кнопок
            btn.setCheckable(True)
            btn.setVisible(False)
            #btn.clicked.connect(self.pr)
            self.hbox_antenna.addWidget(btn)
            btn.setStyleSheet(StyleSheet)
            btn.setObjectName(str(self.list_of_ant_East.index(name)))# устанавливавем имя обьекта равное его индексу в списке
            #btn.clicked.connect(self.onClickbtnAntEast)
            self.lst_btn_ant_East.append(btn) # Формируем список кнопок номеров антенн
    for name in self.list_of_ant_FarEast: # Генерация кнопок по списку антенн
            btn = QtWidgets.QPushButton(name)
            btn.setMinimumWidth(35 )#Устанавливаем минимальную ширину кнопок
            btn.setCheckable(True)
            btn.setVisible(False)
            #btn.clicked.connect(self.pr)
            self.hbox_antenna.addWidget(btn)
            btn.setStyleSheet(StyleSheet)
            btn.setObjectName(str(self.list_of_ant_FarEast.index(name)))# устанавливавем имя обьекта равное его индексу в списке
            #btn.clicked.connect(self.onClickbtnAntFarEast)
            self.lst_btn_ant_FarEast.append(btn) # Формируем список кнопок номеров антенн
                
    self.lst_btn_direction = [] # Создание списка кнопок плечей радиогелиографа
    for direction in self.list_of_directions: # Генерация кнопок по списку плечей радиогелиографа
            btn = QtWidgets.QPushButton(direction)
            btn.setCheckable(True)
            #if direction == 'North':# Добавляем в плечо North центральную антенну
            #    btn.setMinimumWidth(35*9)#Устанавливаем минимальную ширину кнопок
            #else: 
            #    btn.setMinimumWidth(35*8)
            btn.setStyleSheet(StyleSheet)
            self.hbox_direction.addWidget(btn)
            self.lst_btn_direction.append(btn) # Формируем список кнопок номеров антенных групп
    self.lst_btn_direction[0].clicked.connect(self.onClickbtnNorth)
    self.lst_btn_direction[1].clicked.connect(self.onClickbtnWest)
    self.lst_btn_direction[2].clicked.connect(self.onClickbtnEast)
    self.lst_btn_direction[3].clicked.connect(self.onClickbtnFarEast)
    
    #self.show()
    self.move(50,100)
    
    

    
