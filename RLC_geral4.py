import numpy as np
import math
from sympy import symbols, exp, diff, cos, sin
import matplotlib.pyplot as plt

## Resolve o Problema Pratico 8.10
t = symbols('t')

V = float(input('Entre com Entre com a fonte de Tensão: '))
R1 = float(input('Entre com o Resistor R1: '))
R2 = float(input('Entre com o Resistor R2: '))
C1 = float(input('Entre com o Capacitor C1: '))
C2 = float(input('Entre com o Capacitor C2: '))

def linearSol():

	V2 = v2_0
	dv0 = dv2
	coeficientes = [C1*C2*R2**2, C2*R2 + C2*R2 + C1*R1, 1/R1]
	raizes = np.roots(coeficientes)
	s1 = round(raizes[1], 4)
	s2 = round(raizes[0], 4)    
	print('S1:', s1, 'S2', s2)

	if(s1 > s2):
		a = np.array([[1,1],[s1,s2]],dtype='float')
		b = np.array([V2-v2_f,dv0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
		print('A1:',A1, 'A2:',A2)
		resposta = "Amortecimento supercrítico"
	elif(s1 == s2):
		a = np.array([[1,0],[s1, 1]],dtype='float')
		b = np.array([V2-v2_f,dv0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
		print('A1:',A1, 'A2:',A2)
		resposta = "Amortecimento crítico"
	else:
		a = np.array([[1,0],[s1, s2]],dtype='float')
		b = np.array([V2-v2_f,dv0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
		resposta = "Subamortecimento"
	return A1, A2, s1, s2, resposta 


#Para t < 0
v1_0 = 0
v2_0 = 0

print("v1(0) e v2(0):",v1_0,"V")

#Para t = oo
v1_f = V
v2_f = V

print("v1(oo) e v2(oo):",v1_f,"V")

#Para t > 0
dv1 = ((V*R2) - (R2*R1)*v1_0 + (R1*v2_0))/(R1*R2*C1)
dv2 = (v1_0-v2_0)/(R2*C2)

print("dv1(0)/dt:",dv1,"V/s")
print("dv2(0)/dt:",dv2,"V/s")

A1, A2, s1, s2, resposta = linearSol()

if(s1 > s2):
	v2 = v2_f + A1*exp(s1*t) + A2*exp(s2*t)
elif(s1 == s2):
	alpha = s1
	v2 = v2_f + (A1 + A2*t)*exp(-alpha*t)
else:
	alpha = s1.real
	omega_d = s1.imag
	v2 = v2_f + exp(-alpha*t)*(A1*cos(omega_d*t) + A2*sin(omega_d*t))

print("v2(t):",v2, " V")

v1 = C2*diff(v2,t) + v2
print("v1(t):",v1, " V")

vo = v1 - v2
print("########################################")
print("Tipo de Resposta ",resposta)
print("Resposta vo(t):",vo, " V")