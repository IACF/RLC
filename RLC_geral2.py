import numpy as np
import math
from sympy import symbols, exp, diff, cos, sin
import matplotlib.pyplot as plt

t = symbols('t')

circuito = input('Entre com o número do circuito desejato: ')
alpha = 0
omega_d = 0
V = float(input('Entre com a fonte de Tensão: '))
R1 = float(input('Entre com o Resistor R1: '))
R2 = float(input('Entre com o Resistor R2: '))
C = float(input('Entre com o Capacitor C1: '))
L = float(input('Entre com o Indutor: L1 '))

#Para t < 0
i0 = 0
V0 = V

print("i(0):",i0,"A")
print("v(0):",V0,"V")

#Para t = oo
i_f = V/(R1 + R2)
vf = V*R2/(R1 + R2)

print("i(oo):",i_f,"A")
print("v(oo):",vf,"V")

#Para t > 0
ic0 = i0 - (V0/R2)

def linearSol(V0, Vss):

	dv0 = ic0/C

	coeficientes = [L*C, R1*C+(L/R2), 1+(R1/R2)]
	raizes = np.roots(coeficientes)
	s1 = round(raizes[1], 4)
	s2 = round(raizes[0], 4)    
	print('S1:', s1, 'S2', s2)

	if(s1 > s2):
		a = np.array([[1,1],[s1,s2]],dtype='float')
		b = np.array([V0-vf,dv0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
		print('A1:',A1, 'A2:',A2)
		resposta = "Amortecimento supercrítico"
	elif(s1 == s2):
		a = np.array([[1,0],[s1, 1]],dtype='float')
		b = np.array([V0-vf,dv0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
		print('A1:',A1, 'A2:',A2)
		resposta = "Amortecimento crítico"
	else:
		a = np.array([[1,0],[s1, s2]],dtype='float')
		b = np.array([V0-vf,dv0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
		resposta = "Subamortecimento"
	return A1, A2, s1, s2, resposta 


A1, A2, s1, s2, resposta = linearSol(V0, vf)
if(s1 > s2):
	v = A1*exp(s1*t) + A2*exp(s2*t) + vf
elif(s1 == s2):
	v = (A1 + A2*t)*exp(-alpha*t) + vf
else:
	v = exp(-alpha*t)*(A1*cos(omega_d*t) + A2*sin(omega_d*t)) + vf
print("Resposta completa v(t):",v,"V")


i = v/2 + C*diff(v,t)

print("i(t):",i,"A")