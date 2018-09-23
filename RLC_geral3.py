import numpy as np
import math
from sympy import symbols, exp, diff, cos, sin
import matplotlib.pyplot as plt

## Resolve o exemplo 8.10
t = symbols('t')

print("Exemplo 8.10\n")

alpha = 0
omega_d = 0
V = float(input('Entre com a fonte de Tensão: '))
R1 = float(input('Entre com o Resistor R1: '))
R2 = float(input('Entre com o Resistor R2: '))
L1 = float(input('Entre com o Indutor L1: '))
L2 = float(input('Entre com o Indutor L2: '))
I0 = 0

def linearSol():

	If = V/R1
	di0 = V/L1
	coeficientes = [L1*L2, L1 + L2*(R1+R2), R1]
	raizes = np.roots(coeficientes)
	s1 = round(raizes[1], 4)
	s2 = round(raizes[0], 4)    
	print('S1:', s1, 'S2', s2)

	
	# print('A1:',A1, 'A2:',A2)
	if(s1 > s2):
		a = np.array([[1,1],[s1,s2]],dtype='float')
		b = np.array([I0-If,di0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
		resposta = "Amortecimento supercrítico"
	elif(s1 == s2):
		a = np.array([[1,0],[s1, 1]],dtype='float')
		b = np.array([I0-If,di0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
		resposta = "Amortecimento crítico"
	else:
		a = np.array([[1,0],[s1, s2]],dtype='float')
		b = np.array([I0-If,di0],dtype='float')	
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
		resposta = "Subamortecimento"
	return A1, A2, s1, s2, resposta 
#Para t < 0
i1_0 = 0
i2_0 = 0

print("i1(0):",i1_0,"A")
print("i2(0):",i2_0,"A")

#Para t = oo
i_f = V/R1

print("i(oo):",i_f,"A")

#Para t > 0
#di1(0)/dt = vl/L1
di1 = V/L1
#di2(0)/dt = vl/L2
di2 = 0/L2

print("di1(0)/dt:",di1,"A/s")
print("di2(0)/dt:",di2,"A/s")

A1, A2, s1, s2, resposta = linearSol()

print("Raizes s1 e s2: {0} , {1}".format(s1,s2))

print(resposta)
print("Constantes A1 e A2: {0} , {1}".format(A1,A2))

if(s1 > s2):
	i1 = i_f + A1*exp(s1*t) + A2*exp(s2*t)
elif(s1 == s2):
	alpha = s1
	i1 = i_f + (A1 + A2*t)*exp(-alpha*t)
else:
	alpha = s1.real
	omega_d = s1.imag
	i1 = i_f + exp(-alpha*t)*(A1*cos(omega_d*t) + A2*sin(omega_d*t))

print("########################################")
print("Tipo de Resposta ",resposta)
print("i1(t):",i1," A")

i2 = -V + (R1+R2)*i1 + L1*diff(i1,t)
print("i2(t):",i2," A")

vo = i1 - i2

print("V0(t):",vo," V")