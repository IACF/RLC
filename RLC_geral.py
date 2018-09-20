import numpy as np
import math
from sympy import *
import matplotlib.pyplot as plt

t = symbols('t')

circuito = input('Entre com o número do circuito desejato: ')

if(circuito == '1'):

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
elif (circuito == '2'):
	print("Exemplo 8.10\n")

	

	V = 7
	L1 = 1/2
	L2 = 1/5
	R1 = 3
	R2 = 1
	I0 = 0

	def linearSol():

		#Para t = oo
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

	#desligar fontes indep.
	#3i1 + 1/2*di1/dt + (i1 - i2) = 0
		#4i1 + 1/2*di1/dt - i2 = 0
	#1/5*di2/t + i2 - i1 = 0
		#4/5*di1/dt + 1/10*d^2i1/dt^2 + 4i1 + 1/2*di1/dt - i1 = 0
		#d^2i1/dt^2 + 13di1/dt + 30i1 = 0

	#s^2 + 13s + 30 = 0

	# r = solve(s**2 + 13*s + 30,s)
	A1, A2, s1, s2, resposta = linearSol()

	print("Raizes s1 e s2: {0} , {1}".format(s1,s2))

	#raizes reais e negativas: Superamortecido

	#i1(t) = 7/3 + A1*exp(-10t) + A2*exp(-3t)
	#i1(0) = 7/3 + A1 + A2 = 0
		#A1 = -7/3 - A2
	#di1(0)/dt = -10A1 -3A2 = 14
		#-10(-7/3 - A2) - 3A2 = 14
	# A2 = (14 - 70/3)/7
	# A1 = -7/3 - A2
	print(resposta)
	print("Constantes A1 e A2: {0} , {1}".format(A1,A2))

	i1 = i_f + A1*exp(s1*t) + A2*exp(s2*t)

	print("i1(t):",i1,"A")
	
	#V = 3i1 + L1*di1/dt + (i1 - i2)
	# i2 = 3*i1 + L1*diff(i1,t) + i1 - V
	i2 = -V + (R1+R2)*i1 + L1*diff(i1,t)
	print("i2(t):",i2,"A")

	vo = i1 - i2

	print("V0(t):",vo,"V")

elif circuito == '3':
	print("Problema Prático 8.10")

	V = 20
	C1 = 1/2
	C2 = 1/3

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

	#desligar fontes indep.
	#v1/1 + C1*dv1/dt + vo/1 = 0
	#vo = v1-v2
		#v1 + 1/2*dv1/dt + v1-v2 = 0
		#dv1/dt + 4v1 - 2v2 = 0
	#v1 = 1*C2*dv2/dt + v2
		#1/3*d^2v2/dt^2 + dv2/dt + 4/3*dv2/dt + 4v2 - 2v2 = 0
		#d^2v2/dt^2 + 7dv2/dt + 6v2 = 0

	#s^2 + 7s + 6 = 0

	r = solve(s**2 + 7*s + 6,s)
	s1,s2 = r[0],r[1]

	print("Raizes para v2:",s1,s2)

	#raizes reais e negativas: Superamortecido

	#v2(t) = 20 + A1*exp(-6t) + A2*exp(-t)
	#v2(0) = 20 + A1 + A2 = 0
		#A2 = -20 - A1
	#dv2(0)/dt = -6A1 - A2 = 0
		#-6A1 - (-20 - A1) = 0
	A1 = 20/5
	A2 = -20 - A1

	print("Constantes A1 e A2:",A1,A2)

	v2 = v2_f + A1*exp(s1*t) + A2*exp(s2*t)
	print("v2(t):",v2,"V")

	v1 = C2*diff(v2,t) + v2
	print("v1(t):",v1,"V")

	vo = v1 - v2
	print("Resposta vo(t):",vo,"V")