'''
Created on 1 Jul 2018

@author: markot
'''

from PyQt5 import QtWidgets, QtCore
#import logging
import sys
import time
#import threading
#import PyTango as pt
import numpy as np
import pyqtgraph as pq
from scipy.signal import tukey
#from scipy.stats import gennorm
from scipy.special import gamma 
from scipy.stats import moment
from __builtin__ import str
#import scipy as sy

class PlotYYWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent=parent)
        self.setupLayout()
        
    def setupLayout(self):
        
        self.plotYYWidget = pq.PlotWidget(useOpenGL=True)
        self.plot11 = self.plotYYWidget.plot()
        self.plot11.setPen((10, 200, 70))
        self.plot12 = self.plotYYWidget.plot()
        self.plot12.setPen((255,255,0))
        self.plot13 = self.plotYYWidget.plot()
        plot13pen = pq.mkPen(color=(10, 200, 70), width=1, style=QtCore.Qt.DashLine)
        self.plot13.setPen(plot13pen)
        self.plot14 = self.plotYYWidget.plot()
        plot14pen = pq.mkPen(color=(10, 200, 70), width=1, style=QtCore.Qt.DotLine)
        self.plot14.setPen(plot14pen)
        self.plot15 = self.plotYYWidget.plot()
        plot15pen = pq.mkPen(color=(64, 224, 208), width=1, style=QtCore.Qt.DotLine)
        self.plot15.setPen(plot15pen)

        self.plotYYItem = self.plotYYWidget.plotItem
        self.plotYYItem.setLabels(left='abs(El)')
        self.ViewboxYY = pq.ViewBox()
        self.plotYYItem.showAxis('right')
        self.plotYYItem.scene().addItem(self.ViewboxYY)
        self.plotYYItem.getAxis('right').linkToView(self.ViewboxYY)
        self.ViewboxYY.setXLink(self.plotYYItem)
        self.plotYYItem.getAxis('right').setLabel('Phase / rad')
        self.plot21 = pq.PlotCurveItem()
        self.plot21.setPen((200, 70, 10))
        self.ViewboxYY.addItem(self.plot21)
        self.plot22 = pq.PlotCurveItem()
        self.ViewboxYY.addItem(self.plot22)
        plot22pen = pq.mkPen('r', width=1, style=QtCore.Qt.DashLine)
        self.plot22.setPen(plot22pen)
        self.ViewboxYY.addItem(self.plot22)

        #self.plotYYItem.vb.sigResized.connect(self.updateElPlotView)
        self.ViewboxYY.setGeometry(self.plotYYItem.vb.sceneBoundingRect())
        self.ViewboxYY.linkedViewChanged(self.plotYYItem.vb, self.ViewboxYY.XAxis)
        self.plotYYWidget.setAntialiasing(True)
        self.plotYYWidget.showGrid(True, False)
        self.plotYYWidget.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.update()
        
        self.layout = QtWidgets.QHBoxLayout(self)
        self.layout.addWidget(self.plotYYWidget)
        
class SpectrometerCamera(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        #super().__init__()
        
        #self.lock = threading.Lock()       
        #self.scanTimer = QtCore.QTimer()
        #self.scanTimer.timeout.connect(self.updateImage)    
        self.stop_timer = False
        self.secondOrderPhaseCoef = 0.0
        self.thirdOrderPhaseCoef = 0.0
        self.shaperPhaseAmpCoef = 0.0
        self.shaperPhaseWidth = 1.0
        self.shaperPhaseCenter = 262.0
        self.shaperPhaseShapeParam = 2.0
        self.gaussianPhaseCoef = 0.0
        self.ampFilterWidth = 10.0
        self.ampLoss = 1.0
        self.ampMaskOn = False
        self.shaperOn = False
        
        #initial spectrum and fft
        self.wlsize = 1024
        self.shaper_phase = np.zeros(self.wlsize)
        self.wl = np.linspace(266, 258, self.wlsize)
        self.absE_wl_mask = np.empty_like(self.wl)
        self.fftsize = 4096*8
        self.ampMask = np.ones(self.wlsize)
        self.tukey_ind_l = 10
        self.tukey_ind_r = self.wlsize-10
        self.tukey_param = 0.0
        
     
        #self.defineSpectrum()
        #self.spectralToTemporal()
        
        #set up and update widgets
        self.setup_layout()    
        self.updateFFT() 
        #self.generateRandomShaperPhase()
        
    
        
    def gennorm_with_var(self,x,A,beta,alpha,mu):
        try:
            pdf = A*beta/(2*alpha*gamma(1/beta))*np.exp(-(np.abs(x-mu)/alpha)**beta)
            return pdf    
        except ValueError:
            pass    
        
    def defineSpectrum(self):
        # higher to lower wavelengths so the frequencies will be in the right order
        # need for np.interp to work correctly
        print 'defining spectrum'
        self.wl_0 = 262.0
        self.w_wl = self.spectralFWHMSpinbox.value()
        self.absE_wl = np.exp(-4.0*np.log(2)*(self.wl-self.wl_0)**2/self.w_wl**2)
        
        #amp mask
        if self.ampMaskOn:
            #w_tukey = 2.0
            #wl_0 = 262
            #tukey_ind_r = np.max(np.where(self.wl>wl_0-w_tukey/2))
            #tukey_ind_l = np.min(np.where(self.wl<wl_0+w_tukey/2))
            #print tukey_ind_l, tukey_ind_r   
            tukey_lim_r = self.ampMaskLeftSpinBox.value()
            self.tukey_ind_r = np.max(np.where(self.wl>tukey_lim_r))
            tukey_lim_l = self.ampMaskRightSpinBox.value()  
            self.tukey_ind_l = np.min(np.where(self.wl<tukey_lim_l))
            self.tukey_param = self.ampMaskTukeyParamSpinBox.value()               
            self.ampMask = np.hstack((np.zeros(self.tukey_ind_l), tukey(self.tukey_ind_r-self.tukey_ind_l,self.tukey_param), np.zeros(len(self.wl)-self.tukey_ind_r)))
        else:
            self.ampMask = np.ones(self.wlsize)
        
        self.absE_wl_mask = self.absE_wl*self.ampMask
        self.ampLoss = np.sum(self.absE_wl_mask)/np.sum(self.absE_wl)
        #np.multiply(self.absE_wl,self.ampMask,out = self.absE_wl_Mask)
        secondOrderPhase = self.secondOrderPhaseCoef*(self.wl-self.wl_0)**2
        thirdOrderPhase = self.thirdOrderPhaseCoef*(self.wl-self.wl_0)**3
        #self.shaper_phase = self.shaperPhaseAmpCoef*np.exp(-4.0*np.log(2)*(self.wl-self.shaperPhaseCenter)**2/self.shaperPhaseWidth**2)
        self.shaper_phase = self.gennorm_with_var(self.wl,self.shaperPhaseAmpCoef, self.shaperPhaseShapeParam, self.shaperPhaseWidth, self.wl_0)
        if self.shaperOn:
            self.phi_wl = secondOrderPhase + thirdOrderPhase + self.shaper_phase
        else:
            self.phi_wl = secondOrderPhase + thirdOrderPhase          
        print 'done defining spectrum'
                
    def spectralToTemporal(self):
        # lambda -> nu, sample nu uniformly and interpolate
        
        print 'performing FFT'
        self.nu = 3e8/self.wl*1e9/1e12
        nuI = 3e8/self.wl[0]*1e9/1e12
        nuF = 3e8/self.wl[-1]*1e9/1e12
        d_nu = (nuF-nuI)/len(self.wl)
        self.nu_uniform = np.linspace(nuI,nuF,self.wlsize)

        #interpolate
        self.absEnu = np.interp(self.nu_uniform, self.nu, self.absE_wl)
        self.absEnu_mask = np.interp(self.nu_uniform, self.nu, self.absE_wl_mask)
        #self.absEnu = np.multiply(self.absEnu,self.ampMask)
        self.phi_nu=np.interp(self.nu_uniform, self.nu, self.phi_wl)
        self.Enu = np.sqrt(self.absEnu) * np.exp(1j*self.phi_nu)
        self.Enu_mask = np.sqrt(self.absEnu_mask) * np.exp(1j*self.phi_nu)

        pad_length = np.round((self.fftsize - len(self.nu_uniform))/2)
        self.Enu_padded = np.pad(self.Enu,(pad_length,pad_length),'constant')
        self.Enu_mask_padded = np.pad(self.Enu_mask,(pad_length,pad_length),'constant')

        self.Et = np.fft.fftshift(np.fft.fft(np.fft.fftshift(self.Enu_padded)))
        self.Et_mask = np.fft.fftshift(np.fft.fft(np.fft.fftshift(self.Enu_mask_padded)))
        #Et=np.fft.fftshift(np.fft.fft(abs(self.Enu_padded)))
        #Et=np.fft.fft(abs(self.Enu_padded))
        self.t=np.linspace(-1/2.0/d_nu,1/2.0/d_nu,self.fftsize)
        self.Intt=np.abs(self.Et)**2/np.max(abs(self.Et)**2)
        self.Intt_mask = np.abs(self.Et_mask)**2/np.max(abs(self.Et_mask)**2)
        self.phit=np.unwrap(np.angle(self.Et))
        self.phit_mask = np.unwrap(np.angle(self.Et_mask))
        print 'done with FFT'
    
    def arrayThreshold(self,a):
        #print 'thresholding array'
        indices = np.where(a > np.max(a)/10000.0)
        #indices = np.where(a > 0.0)
        ind_left = np.min(indices)
        ind_right = np.max(indices)
        return ind_left, ind_right
    
    def updateElPlotView(self):
        self.ViewboxYY.setGeometry(self.plotYYItem.vb.sceneBoundingRect())
        self.ViewboxYY.linkedViewChanged(self.plotYYItem.vb, self.ViewboxYY.XAxis)                
        
    def updateFFT(self):
        print 'updating FFT plots'
        self.secondOrderPhaseCoef = self.secondOrderPhaseCoefSpinbox.value()
        self.thirdOrderPhaseCoef = self.thirdOrderPhaseCoefSpinbox.value()
        self.shaperPhaseAmpCoef = self.shaperPhaseAmpSpinbox.value()
        self.shaperPhaseWidth = self.shaperPhaseWidthSpinbox.value()
        self.shaperPhaseShapeParam = self.shaperPhaseShapeParamSpinbox.value()
        #print self.secondOrderPhaseCoef
        
        self.defineSpectrum()
        self.spectralToTemporal()
        
        self.plotSpectralWidget.plot11.setData(self.wl,self.absE_wl)
        self.plotSpectralWidget.plot12.setData(self.wl,self.ampMask)
        self.plotSpectralWidget.plot21.setData(self.wl,self.phi_wl)
        self.plotSpectralWidget.plot22.setData(self.wl,self.shaper_phase)
        
        [self.ind1, self.ind2] = self.arrayThreshold(self.Intt)
        self.plotTemporalWidget.plot11.setData(self.t[self.ind1:self.ind2],self.Intt[self.ind1:self.ind2])
        self.plotTemporalWidget.plot13.setData(self.t[self.ind1:self.ind2],self.Intt_mask[self.ind1:self.ind2])
        self.plotTemporalWidget.plot21.setData(self.t[self.ind1:self.ind2],self.phit[self.ind1:self.ind2]-np.mean(self.phit[self.ind1:self.ind2])) 
        
        print 'done updating FFT plots'
        #if self.stop_timer is not True:
        #    self.cameraTimer.start(100)
            
    def closeEvent(self, event):
        self.stop_timer = True
        self.deleteLater()
        time.sleep(0.1)
        print 'stopping'
        
    def updateEverything(self):
        print 'updating everything'
        self.defineSpectrum()
        self.spectralToTemporal()
        self.updateFFT()
        #print str(self.ampLoss,'%6.2f')
        
    def toggleAmpMask(self):
        if self.ampFilterRadioButton.isChecked():
            self.ampMaskOn = True
        else:
            self.ampMaskOn = False
        self.updateEverything()
        
    def applyShaperPhase(self):
        print 'shaper phase toggled'
        if self.shaperRadioButton.isChecked():
            self.shaperOn = True
        else:
            self.shaperOn = False
        self.updateEverything()    
                    
    def generateRandomShaperPhase(self):
        noGenes = 25
        genAmp=2
        genes = (np.random.rand(noGenes)*2-1)*genAmp
        genes -= genes[0]
        phase_coarse = np.cumsum(genes)-genes            
        phase = np.interp(np.linspace(0,noGenes-1,len(self.nu_uniform)), range(noGenes), phase_coarse)
        #print phase_coarse     
        self.shaper_phase = phase
        self.updateEverything()
        
    def exportPulseShape(self):
        fname = self.exportFileName.text()
        print fname
        spectral = np.array([self.wl, self.absE_wl, self.ampMask, self.phi_wl, self.shaper_phase]);
        temporal = np.array([self.t[self.ind1:self.ind2],self.Intt[self.ind1:self.ind2],self.Intt_mask[self.ind1:self.ind2],self.phit[self.ind1:self.ind2]-np.mean(self.phit)])
        np.savetxt(fname + 'spectral.txt', spectral)
        np.savetxt(fname + 'temporal.txt', temporal)
    
    def setup_layout(self):
        
        print 'setting up layout'
        self.layout = QtWidgets.QHBoxLayout(self) #the whole window, main layout
        self.plotLayout = QtWidgets.QVBoxLayout()
        self.controlsLayout = QtWidgets.QVBoxLayout()

        self.plotSpectralWidget = PlotYYWidget()
        #self.plotSpectralWidget.plot11.setData(self.wl,self.absE_wl)
        #self.plotSpectralWidget.plot12.setData(self.wl,self.ampMask)
        #self.plotSpectralWidget.plot21.setData(self.wl,self.phi_wl)
        #self.plotSpectralWidget.plot22.setData(self.wl,self.shaper_phase)
        self.plotSpectralWidget.show()
        self.plotLayout.addWidget(self.plotSpectralWidget)
        
        self.plotTemporalWidget = PlotYYWidget()
        #[ind1, ind2] = self.arrayThreshold(self.Intt)
        #self.plotTemporalWidget.plot11.setData(self.t[ind1:ind2],self.Intt[ind1:ind2])
        #self.plotTemporalWidget.plot13.setData(self.t[ind1:ind2],self.Intt_mask[ind1:ind2])
        #self.plotTemporalWidget.plot21.setData(self.t[ind1:ind2],self.phit[ind1:ind2]-np.mean(self.phit[ind1:ind2]))
        self.plotLayout.addWidget(self.plotTemporalWidget)       
           
#         self.spectralPhaseSlider = mySlider()
#         self.spectralPhaseSlider.slider.setValue(self.secondOrderPhaseCoef)
#         self.spectralPhaseSlider.label.setText(str(self.secondOrderPhaseCoef))
#         self.spectralPhaseSlider.valueChanged.connect(self.updateFFTplots)
#         self.controlsLayout.addWidget(self.spectralPhaseSlider)

        self.spectralFWHMSpinbox = QtWidgets.QDoubleSpinBox()
        self.spectralFWHMSpinbox.setDecimals(2)
        self.spectralFWHMSpinbox.setValue(1.6)
        self.spectralFWHMSpinbox.setMinimum(0.1)
        self.spectralFWHMSpinbox.valueChanged.connect(self.updateEverything)
        self.spectralFWHMLabel = QtWidgets.QLabel('spectral width')

        self.secondOrderPhaseCoefSpinbox = QtWidgets.QDoubleSpinBox()
        self.secondOrderPhaseCoefSpinbox.setDecimals(2)
        self.secondOrderPhaseCoefSpinbox.setValue(self.secondOrderPhaseCoef)
        self.secondOrderPhaseCoefSpinbox.valueChanged.connect(self.updateEverything)
        self.secondOrderPhaseCoefSpinbox.setMinimum(-150)
        self.secondOrderPhaseCoefSpinbox.setMaximum(150)
        self.secondOrderPhaseCoefLabel = QtWidgets.QLabel()
        self.secondOrderPhaseCoefLabel.setText('2nd order phase coef')
        vspaceritem = QtWidgets.QSpacerItem(100,100,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding)           
           
        self.thirdOrderPhaseCoefSpinbox = QtWidgets.QDoubleSpinBox()
        self.thirdOrderPhaseCoefSpinbox.setDecimals(2)
        self.thirdOrderPhaseCoefSpinbox.setValue(self.thirdOrderPhaseCoef)
        self.thirdOrderPhaseCoefSpinbox.valueChanged.connect(self.updateEverything)
        self.thirdOrderPhaseCoefSpinbox.setMinimum(-150)
        self.thirdOrderPhaseCoefSpinbox.setMaximum(150)
        self.thirdOrderPhaseCoefLabel = QtWidgets.QLabel()
        self.thirdOrderPhaseCoefLabel.setText('3rd order phase coef')           
           
        self.controlsLayout.addWidget(self.spectralFWHMLabel)
        self.controlsLayout.addWidget(self.spectralFWHMSpinbox)    
        self.controlsLayout.addWidget(self.secondOrderPhaseCoefLabel)
        self.controlsLayout.addWidget(self.secondOrderPhaseCoefSpinbox)
        self.controlsLayout.addWidget(self.thirdOrderPhaseCoefLabel)
        self.controlsLayout.addWidget(self.thirdOrderPhaseCoefSpinbox)
        self.controlsLayout.addItem(vspaceritem)
        
        # shaper controls
        self.shaperRadioButton = QtWidgets.QCheckBox()
        self.shaperRadioButton.toggled.connect(self.applyShaperPhase)
        self.shaperRadioButtonLabel = QtWidgets.QLabel('apply shaper phase')
        
        self.shaperPhaseAmpSpinbox = QtWidgets.QDoubleSpinBox()
        self.shaperPhaseAmpSpinbox.setDecimals(2)
        self.shaperPhaseAmpSpinbox.setValue(self.shaperPhaseAmpCoef)
        self.shaperPhaseAmpSpinbox.valueChanged.connect(self.updateEverything)
        self.shaperPhaseAmpSpinbox.setMinimum(-50)
        self.shaperPhaseAmpSpinbox.setMaximum(50)
        self.shaperPhaseAmpLabel = QtWidgets.QLabel()
        self.shaperPhaseAmpLabel.setText('shaper phase amp.')

        self.shaperPhaseShapeParamSpinbox = QtWidgets.QDoubleSpinBox()
        self.shaperPhaseShapeParamSpinbox.setDecimals(2)
        self.shaperPhaseShapeParamSpinbox.setValue(self.shaperPhaseShapeParam)
        self.shaperPhaseShapeParamSpinbox.editingFinished.connect(self.updateEverything)
        self.shaperPhaseShapeParamLabel = QtWidgets.QLabel()
        self.shaperPhaseShapeParamLabel.setText('shaper phase shape param.')       
        
        self.shaperPhaseWidthSpinbox = QtWidgets.QDoubleSpinBox()
        self.shaperPhaseWidthSpinbox.setDecimals(2)
        self.shaperPhaseWidthSpinbox.setValue(self.shaperPhaseWidth)
        self.shaperPhaseWidthSpinbox.editingFinished.connect(self.updateEverything)
        self.shaperPhaseWidthSpinbox.setMinimum(0)
        self.shaperPhaseWidthSpinbox.setMaximum(6)
        self.shaperPhaseWidthLabel = QtWidgets.QLabel()
        self.shaperPhaseWidthLabel.setText('shaper phase width')
        
        #self.applyRandomPhaseButton = QtWidgets.QPushButton('Apply random phase')
        #self.applyRandomPhaseButton.clicked.connect(self.generateRandomShaperPhase)
        
        self.controlsLayout.addWidget(self.shaperRadioButtonLabel)
        self.controlsLayout.addWidget(self.shaperRadioButton)
        self.controlsLayout.addWidget(self.shaperPhaseAmpLabel)
        self.controlsLayout.addWidget(self.shaperPhaseAmpSpinbox)
        self.controlsLayout.addWidget(self.shaperPhaseShapeParamLabel)
        self.controlsLayout.addWidget(self.shaperPhaseShapeParamSpinbox)
        self.controlsLayout.addWidget(self.shaperPhaseWidthLabel)
        self.controlsLayout.addWidget(self.shaperPhaseWidthSpinbox)   
        #self.controlsLayout.addWidget(self.applyRandomPhaseButton)    

        self.controlsLayout.addItem(vspaceritem)
        
        #amp mask controls
        self.ampFilterRadioButton = QtWidgets.QRadioButton()
        self.ampFilterRadioButtonLabel = QtWidgets.QLabel()
        self.ampFilterRadioButtonLabel.setText('apply amp mask?')
        self.ampMaskWidthlabel = QtWidgets.QLabel('amp. mask width')
        self.ampMaskWidthSpinBox = QtWidgets.QDoubleSpinBox()
        self.ampMaskWidthSpinBox.valueChanged.connect(self.updateEverything)
        self.ampFilterRadioButton.toggled.connect(self.toggleAmpMask)
        self.ampMaskLeftSpinBox = QtWidgets.QDoubleSpinBox()
        self.ampMaskLeftSpinBox.setMinimum(250)
        self.ampMaskLeftSpinBox.setMaximum(300)
        self.ampMaskLeftSpinBox.setValue(np.min(self.wl)+1)
        self.ampMaskLeftSpinBox.valueChanged.connect(self.updateEverything)
        self.ampMaskLeftLabel = QtWidgets.QLabel('mask left limit')        
        self.ampMaskRightSpinBox = QtWidgets.QDoubleSpinBox()
        self.ampMaskRightSpinBox.setMinimum(250)
        self.ampMaskRightSpinBox.setMaximum(300)
        self.ampMaskRightSpinBox.setValue(np.max(self.wl)-1)
        self.ampMaskRightSpinBox.valueChanged.connect(self.updateEverything)
        self.ampMaskRightLabel = QtWidgets.QLabel('mask right limit')
        self.ampMaskTukeyParamSpinBox = QtWidgets.QDoubleSpinBox()
        self.ampMaskTukeyParamSpinBox.setMinimum(0)
        self.ampMaskTukeyParamSpinBox.setMaximum(1)
        self.ampMaskTukeyParamSpinBox.setValue(0)
        self.ampMaskTukeyParamSpinBox.valueChanged.connect(self.updateEverything)
        self.ampMaskTukeyParamLabel = QtWidgets.QLabel('Tukey param.')
        self.ampMaskLossLabel = QtWidgets.QLabel()
        self.ampMaskLossLabel.setText(str(self.ampLoss))
        self.controlsLayout.addWidget(self.ampFilterRadioButtonLabel)
        self.controlsLayout.addWidget(self.ampFilterRadioButton)
        #self.controlsLayout.addWidget(self.ampMaskWidthlabel)
        #self.controlsLayout.addWidget(self.ampMaskWidthSpinBox)
        self.controlsLayout.addWidget(self.ampMaskLeftLabel)
        self.controlsLayout.addWidget(self.ampMaskLeftSpinBox)
        self.controlsLayout.addWidget(self.ampMaskRightLabel)
        self.controlsLayout.addWidget(self.ampMaskRightSpinBox)
        self.controlsLayout.addWidget(self.ampMaskTukeyParamLabel)
        self.controlsLayout.addWidget(self.ampMaskTukeyParamSpinBox)  
        self.controlsLayout.addWidget(self.ampMaskLossLabel)      
        self.controlsLayout.addItem(vspaceritem)
        

        
        self.exportFileName = QtWidgets.QLineEdit()
        self.exportButton = QtWidgets.QPushButton('Export')
        self.exportButton.clicked.connect(self.exportPulseShape)
        self.controlsLayout.addWidget(self.exportFileName)
        self.controlsLayout.addWidget(self.exportButton)
                
        self.layout.addLayout(self.plotLayout)    
        self.layout.addLayout(self.controlsLayout)   
                  
#         self.plotWidget3 = pq.PlotWidget(useOpenGL=True)                
#         self.plotWidget3.setSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Maximum)
#         self.plot31 = self.plotWidget3.plot()
#         self.plot31.setData(self.arrayThreshold(self.Intt, self.t)[1],self.arrayThreshold(self.Intt, self.t)[0])
#         self.plot32 = self.plotWidget3.plot()
#         self.plot32.setData(self.t,self.Intt,color='r')
#         self.plot32.setPen('r')
        #self.layout.addWidget(self.plotWidget3) 
        
class mySlider(QtWidgets.QSlider, QtWidgets.QLabel):
    #doesn't work yet!!!
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent=parent)
        self.setupLayout()
        
    def setupLayout(self):
        self.slider = QtWidgets.QSlider()
        self.slider.setGeometry(10, 10, 10, 10)
        self.slider.setTickInterval(1)
        self.slider.setMaximum(50)
        self.slider.setMinimum(-50)
        #self.slider.setMinimumWidth(40)
        self.slider.setSizePolicy(QtWidgets.QSizePolicy.Minimum,QtWidgets.QSizePolicy.Minimum)
        self.slider.setTickPosition(QtWidgets.QSlider.TicksLeft)
        
        self.label = QtWidgets.QLabel()
        self.label.setText('test label')
        
        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.slider)
        self.layout.addWidget(self.label)       
        
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    #app.processEvents()
    myapp = SpectrometerCamera()
    myapp.show()
    #splash.finish(myapp)
    sys.exit(app.exec_())
