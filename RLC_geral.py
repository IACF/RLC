import numpy as np
import math
from sympy import *
import matplotlib.pyplot as plt

t = symbols('t')

V = 12
C = 1/2
L = 1
R1 = 4
R2 = 2

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
#desativar fontes independentes

#i = v/R2 + C*dv/dt

#4i + L*di/dt + v = 0
	#4*(v/2 + 1/2*dv/dt) + d(v/2 + 1/2*dv/dt)/dt + v = 0
	#2v + 2dv/dt + 1/2*dv/dt + 1/2*d^2v/t^2 + v = 0
	#d^2v/dt^2 + 5dv/dt + 6v = 0
	
#s^2 + 5s + 6 = 0

def linearSol(V0, Vss):

	dv0 = ic0/C

	coeficientes = [L*C, R1*C+(L/R2), 1+(R1/R2)]
	raizes = np.roots(coeficientes)
	s1 = round(raizes[1], 4)
	s2 = round(raizes[0], 4)    
	print('S1:', s1, 'S2', s2)

	a = np.array([[1,1],[s1,s2]],dtype='float')
	b = np.array([V0-vf,dv0],dtype='float')
	x = np.linalg.solve(a,b)
	A1 = x[0]
	A2 = x[1]
	print('A1:',A1, 'A2:',A2)
	if(s1 > s2):
		resposta = "Amortecimento supercrítico"
	elif(s1 == s2):
		resposta = "Amortecimento crítico"
	else:
		resposta = "Subamortecimento"
	return A1, A2, s1, s2, resposta 



#Raizes reais e negativas: Superamortecido

#vn(t) = A1*exp(-2t) + A2*exp(-3t)
#vss(t) = v(oo) = 4

#v(t) = 4 + A1*exp(-2t) + A2*exp(-3t)

#dv(0)/dt = -2A1 -3A2 = ic(0)/C
	#ic(0) = -6
	#C = 1/2
	#-2A1 - 3A2 = -12
	#2A1 + 3A2 = 12

#v(0) = 4 + A1 + A2 = 12
	#A1 + A2 = 8

#2(8 - A2) + 3A2 = 12
# A2 = -4.
# A1 = 12.

A1, A2, s1, s2, resposta = linearSol(V0, vf)
v = A1*exp(s1*t) + A2*exp(s2*t) + vf

print("Resposta completa v(t):",v,"V")

#i = v/2 + C*dv/dt

i = v/2 + C*diff(v,t)

print("i(t):",i,"A")