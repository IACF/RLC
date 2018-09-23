import numpy as np
import math
from sympy import symbols, exp, diff
import matplotlib.pyplot as plt

t = symbols('t')

circuito = input('Entre com o número do circuito desejato: ')

if(circuito == '1'):
	print("Exemplo 8.9\n")


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


	A1, A2, s1, s2, resposta = linearSol(V0, vf)
	v = A1*exp(s1*t) + A2*exp(s2*t) + vf

	print("Resposta completa v(t):",v,"V")


	i = v/R2 + C*diff(v,t)

	print("i(t):",i,"A")


elif (circuito == '2'):
	print("Exemplo 8.10\n")

	V = float(input('Entre com a fonte de Tensão: '))
	R1 = float(input('Entre com o Resistor R1: '))
	R2 = float(input('Entre com o Resistor R2: '))
	L1 = float(input('Entre com o Indutor L1: '))
	L2 = float(input('Entre com o Indutor L2: '))
	I0 = 0

	def linearSol():

		If = V/R1
		di0 = V/L1
		print('A:',L1*L2, 'B:',L1+L2*(R1*R2) ,'C', R1)
		coeficientes = [L1*L2, L1 + L2*(R1+R2), R1]
		raizes = np.roots(coeficientes)
		s1 = round(raizes[1], 4)
		s2 = round(raizes[0], 4)    
		print('S1:', s1, 'S2', s2)

		a = np.array([[1,1],[s1,s2]],dtype='float')
		b = np.array([I0-If,di0],dtype='float')
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

	i1 = i_f + A1*exp(s1*t) + A2*exp(s2*t)

	print("i1(t):",i1,"A")
	
	i2 = -V + (R1+R2)*i1 + L1*diff(i1,t)
	print("i2(t):",i2,"A")

	vo = i1 - i2

	print("V0(t):",vo,"V")

elif circuito == '3':
	print("Problema Prático 8.10")


	V = float(input('Entre com Entre com a fonte de Tensão: '))
	R1 = float(input('Entre com o Resistor R1: '))
	R2 = float(input('Entre com o Resistor R2: '))
	C1 = float(input('Entre com o Capacitor C1: '))
	C2 = float(input('Entre com o Capacitor C2: '))
	L = float(input('Entre com o Indutor L1: '))

	def linearSol():

		V2 = 0
		If = V/R1
		dv0 = 0
		coeficientes = [C1*C2*R2**2, C2*R2 + C2*R2 + C1*R1, 1/R1]
		raizes = np.roots(coeficientes)
		s1 = round(raizes[1], 4)
		s2 = round(raizes[0], 4)    
		print('S1:', s1, 'S2', s2)

		a = np.array([[1,1],[s1,s2]],dtype='float')
		b = np.array([V2-V,dv0],dtype='float')
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


	#Para t < 0
	v1_0 = 0
	v2_0 = 0

	print("v1(0) e v2(0):",v1_0,"V")

	#Para t = oo
	v1_f = V
	v2_f = V

	print("v1(oo) e v2(oo):",v1_f,"V")

	#Para t > 0
	#dv1(0)/dt = i1(0)/C1 = (V/1)/(1/2)
	dv1 = V/C1
	#dv2(0)/dt = i2(0)/C2 = 0/C2
	dv2 = 0

	print("dv1(0)/dt:",dv1,"V/s")
	print("dv2(0)/dt:",dv2,"V/s")

	A1, A2, s1, s2, resposta = linearSol()

	print("Constantes A1 e A2:",A1,A2)

	v2 = v2_f + A1*exp(s1*t) + A2*exp(s2*t)
	print("v2(t):",v2,"V")

	v1 = C2*diff(v2,t) + v2
	print("v1(t):",v1,"V")

	vo = v1 - v2
	print("Resposta vo(t):",vo,"V")

print('\n')
input('Digite enter para finalizar .....')
