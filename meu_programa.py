from PyQt5 import QtCore, QtGui, QtWidgets
import sys
import numpy as np
import math
from sympy import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from RLC_novo import linearSol, resposta_rlc, imprime_resultado
import inter
import random

t = symbols('t')
 
class Main(QtWidgets.QMainWindow, inter.Ui_MainWindow):
    def __init__(self, parent=None):
        super(Main, self).__init__(parent)
        self.setupUi(self)

        self.figure = plt.figure()

        self.canvasGraph = FigureCanvas(self.figure)

        # self.toolbar = NavigationToolbar(self.canvas, self)

        #self.connect(self.gerarResultados, SIGNAL("clicked()"), self.imprimeResultados)
        self.gerarResultados.clicked.connect(self.imprimeResultados)
        self.gerarResultados.clicked.connect(self.plot)

    def plot(self):
        data = [random.random() for i in range(10)]

        self.figure.clear()

        ax = self.figure.add_subplot(111)

        ax.plot(data, '*-')

        # refresh canvas
        self.canvasGraph.draw()

    def imprimeResultados(self):
        R = float(self.lineEdit_resistencia.text())
        C = float(self.lineEdit_capacitor.text())
        L = float(self.lineEdit_Indutor.text())
        V0 = float(self.lineEdit_V0.text())
        I0 = float(self.lineEdit_I0.text())

        #print(self.comboBox.currentText())
        #print(self.comboBox.currentIndex())
        #print(self.comboBox.findText('RLC em Série sem Fonte'))

        if(self.comboBox.currentText() == 'RLC em Série sem Fonte'):
            print('RLC em serie')
            associacao = '--'
            alpha = float(R/(2*L))
            omega = float(1./(sqrt(L*C)))
            A1, A2, s1, s2 = linearSol(alpha,omega,V0, I0, associacao, R, C, L)
            resposta,r,r_degrau = resposta_rlc(alpha, omega, associacao, s1, s2, A1, A2, 0, 0)
            print(resposta)
            ylabel = imprimeTodosResultados(self, alpha, omega, resposta, 0, r, r_degrau, I0, V0, L, associacao)
        
            #plotaGrafico(ylabel)
            #self.message.clear()
            #saida = "Tipo de Resposta ", resposta +'\n' + str(r)
            #self.message.setPlainText("saidaaaaaaaaaaaaaaaa")
            #self.message.insertPlainText("Tipo de Resposta " + resposta )
            # self.message.appendPlainText("Tipo de Resposta " + resposta + '\n')
            #self.message.appendPlainText("Resposta i(t): " + str(r) + " A")
            imprime_resultado(r, resposta, alpha, omega, I0, V0, L, 0, r_degrau, associacao) 
    

        


        print('botao clicado ', self.lineEdit_resistencia.text())
        #self.message.setPlainText(self.lineEdit_capacitor.text())

# def inicializaCanvas(self, img):

#     tx = np.arange(0.,3.5,0.1)
#     f = lambdify(t,r)  # converte expressão em função
#     ty = f(tx)
#     plt.xlabel('tempo (s)')
# 	#plt.ylabel(ylabel)
#     #self.plotWidget.canvas.ax.plot()
#     # fig = Figure()
#     # fig.clear()
#     # Graph = fig.

def imprimeTodosResultados(self, alpha, omega, resposta, degrau, r, r_degrau, I0, V0, L, associacao):
    self.message.clear()
    self.message.insertPlainText("Alpha: " + str("{0:10.5f}".format(alpha)) + " Np/s" + '\n' + 
                                "Omega: " + str("{0:10.5f}".format(omega)) + " rad/s")
    self.message.appendPlainText("########################################")
    self.message.appendPlainText("Tipo de Resposta: " + resposta)
    if(I0 == 0):
        self.message.appendPlainText("Resposta i(t): " + str(r) + " A")
        ylabel = 'corrente (A)'				# Titulo da Ordenada do grafico
        if(degrau == '1'):
            self.message.appendPlainText("Resposta ao degrau i(t): "  + str(r_degrau) + " A")
        self.message.appendPlainText("########################################")
    elif(V0 == 0):
        if(alpha == 0):   # restriçao para um caso em especifico
            r = r.diff(t)*(-L)
        self.message.appendPlainText("Resposta v(t): " + str(r) + " V")
        ylabel = 'tensão (V)'
        if(degrau == '1'):
            self.message.appendPlainText("Resposta ao degrau v(t): "  + str(r_degrau) + " V")
        self.message.appendPlainText("########################################")
    elif(associacao == '--'):
        self.message.appendPlainText("Resposta i(t): " + str(r) + " A")
        ylabel = 'corrente (A)'				# Titulo da Ordenada do grafico
        if(degrau == '1'):
            self.message.appendPlainText("Resposta ao degrau i(t): "  + str(r_degrau) + " A")
        self.message.appendPlainText("########################################")
    elif(associacao == '||'):
        self.message.appendPlainText("Resposta v(t): " + str(r) + " V")
        ylabel = 'tensão (V)'
        if(degrau == '1'):
            self.message.appendPlainText("Resposta ao degrau v(t): "  + str(r_degrau) + " V")
        self.message.appendPlainText("########################################")
    return ylabel    

def main():    
    app = QtWidgets.QApplication(sys.argv)
    main_window = Main()
    main_window.show()
    app.exec_()
 
if __name__ == '__main__':
    main()