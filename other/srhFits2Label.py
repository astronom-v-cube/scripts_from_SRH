# -*- coding: utf-8 -*-
"""
Created on Mon May 30 02:54:01 2016

@author: Sergey
"""
import sys, numpy
from PyQt5 import QtGui
from PyQt5 import QtWidgets
import srhFitsFile


class QLabelWithMouse(QtWidgets.QLabel):
    def __init__(self,parent=None):
        QtWidgets.QLabel.__init__(self,parent)

    def mouseMoveEvent(self,event):
        xx = event.x()
        yy = event.y()
        status.setText('X %d, Y %d, EW %d, S %d'%(xx,yy,xx/8 + 49,abs(yy/8) + 177))
        
app = QtWidgets.QApplication.instance();
if not app:
    app = QtWidgets.QApplication(sys.argv);

sF = srhFitsFile.SrhFitsFile();
sF.setCalibIndex(10);
sF.open('/home/sergey/SRH/2017/05/29/mf_20170529_030201.fit');
sF.vis2uv(2,30,phaseCorrect=True);
sF.uv2lmImage(2,30);
#data = numpy.abs(sF.uvRcp)
data = sF.lcp * .1
#data = numpy.clip(data,0.,0.1*data.max())
im = (data - data.min())*255./(data.max() - data.min())
#im = im[127 - 32:127 + 32,127 - 32:127 + 32]
#im = numpy.repeat(numpy.repeat(im,4,axis=0),4,axis=1)
im = numpy.require(im, numpy.uint8, 'C')

image = QtGui.QImage(im.data,im.shape[1],im.shape[0],im.strides[0],QtGui.QImage.Format_Indexed8)

canv = QLabelWithMouse(app.activeWindow())
#canv.setPixmap(QtGui.QPixmap('i160531_0315.png'))
canv.setPixmap(QtGui.QPixmap(image))

status = QtWidgets.QLabel(app.activeWindow())

centralWidget = QtWidgets.QWidget(app.activeWindow());
#centralWidget.setGeometry(0,15,300,290);

layout = QtWidgets.QVBoxLayout(centralWidget);
layout.addWidget(canv);
layout.addWidget(status);     

#canv.show()
#status.show()

centralWidget.show()
sys.exit(app.exec_())
