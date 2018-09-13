import ahkab
import numpy as np
import math
from sympy import *
import matplotlib.pyplot as plt

# R = float(input('Entre com o Resistor: '))
# C = float(input('Entre com o Capacitor: '))
# V0 = float(input('Entre com a tensão V(0) do Capacitor: '))
# I0 = float(input('Entre com a corrente I(0) do Indutor: '))
I0 = 0
V0 = 5
# Vs = 24
# L = float(input('Entre com o Indutor: '))
associacao = input('Entre com  \"\\\\" para RLC paralelo ou \"--\" para RLC série: ')

m = 10**(-3) ##definicao de mili
R = 6.25
C = 10*m
print('CAPAC:',C)
L = 1

A1 = symbols('A1')
A2 = symbols('A2')
B1 = symbols('B1')
B2 = symbols('B2')
t = symbols('t')

def linearSol(alpha, omega, V0, I0, associacao):

	di0 = -(1/L)*(R*I0+V0)
	dv0 = -(R*I0 + V0)/(R*C)
	print('di0:', dv0)
	if(omega > alpha):
		s1 = -alpha
		s2 = sqrt(abs(alpha**2 - omega**2))
		a = np.array([[1,0],[s1, s2]],dtype='float')
		if(associacao == '--'):
			b = np.array([I0,di0],dtype='float')
		else:
			b = np.array([V0,dv0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
	else:
		s1 = -alpha + sqrt(alpha**2 - omega**2)
		s2 = -alpha - sqrt(alpha**2 - omega**2)
		a = np.array([[1,1],[s1, s2]],dtype='float')
		b = np.array([V0-Vs,dv0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[1]
		A2 = x[0]
	print('A1:',A1, 'A2:',A2)
	return A1, A2, s1, s2





def resposta_rlc(alpha, omega, associacao, Vs, Is, s1, s2, A1, A2): #funcao para verificar tipo de resposta e retornar a resposta natural a partir do tipo de amortecimento
	resposta = ""
	
	if alpha > omega:
		resposta = "supercrítico"
		if(associacao == '\\\\'):
			r = A1*exp(s1*t) + A2*exp(s2*t)
			r_degrau = Is + A1*exp(s1*t) + A2*exp(s2*t)
		elif (associacao == '--') :
			r = A1*exp(s1*t) + A2*exp(s2*t)
			r_degrau = Vs + A1*exp(s1*t) + A2*exp(s2*t)
	elif alpha == omega:
		resposta = "amortecimento critico"
		if(associacao == '\\\\'):
			r = (A1 + A2*t)*exp(-alpha*t)
			r_degrau = Is + (A1 + A2*t)*exp(-alpha*t)
		elif (associacao == '--') :
			r = (A2 + A1*t)*exp(-alpha*t)
			r_degrau = Vs + (A2 + A1*t)*exp(-alpha*t)
	else:
		resposta = "subamortecimento"
		omega_d = sqrt(abs(omega**2 - alpha**2))
		if(associacao == '\\\\'):
			r = exp(-alpha*t)*(A1*cos(omega_d*t) + A2*sin(omega_d*t))
			r_degrau = Is + exp(-alpha*t)*(A1*cos(omega_d*t) + A2*sin(omega_d*t))
		elif (associacao == '--') :
			r = exp(-alpha*t)*(A1*cos(omega_d*t) + A2*sin(omega_d*t))
			r_degrau = exp(-alpha*t)*(A1*cos(omega_d*t) + A2*sin(omega_d*t))
	return resposta,r,r_degrau

#----------------------------------------------------------------------------------------#
if associacao == '--':
	while(True):
		try:
			Vs = float (input('Entre com a tensão Vss do capacitor: '))
			break
		except:
			print('Entrada inválida')


	alpha = R/(2*L)
	omega = 1./(sqrt(L*C))
	# if(omega > alpha):
	# 	s1 = -alpha
	# 	s2 = sqrt(abs(alpha**2 - omega**2))
	# else:
	# 	s1 = -alpha + sqrt(alpha**2 - omega**2)
	# 	s2 = -alpha - sqrt(alpha**2 - omega**2)

	
	# di0 = -(1/L)*(R*I0+V0)
	# dv0 = 16
	# print('di0:', di0)
	# a = np.array([[1,1],[s1, s2]],dtype='float')
	# b = np.array([V0-Vs,dv0],dtype='float')
	# x = np.linalg.solve(a,b)
	# A1 = x[1]
	# A2 = x[0]
	A1, A2, s1, s2 = linearSol(alpha,omega,V0, I0, associacao)


	resposta,r,r_degrau = resposta_rlc(alpha, omega,associacao, Vs, 0, s1, s2, A1, A2)

	print("Alpha:",alpha)
	print("Omega:",omega)
	print("Raiz s1:",s1)
	print("Raiz s2:",s2)
	print("Resposta:",resposta)
	print("Resposta i(t):",r)
	print("Resposta ao degrau i(t):",r_degrau)
	tx = np.arange(0,5,0.1)
	f = lambdify(t,r_degrau)
	ty = f(tx)
	plt.plot(tx,ty)
	plt.show()
	print(f(1))

elif associacao == '\\\\':
	
	while(True):
		try:
			Is = float (input('Entre com I(0) para obter a resposta ao degrau, ou 0 caso não queira: '))
			break
		except:
			print('Entrada inválida')
	alpha = 1/(2*R*C)
	omega = 1/(np.sqrt(L*C))
	# s1 = -alpha + sqrt(alpha**2 - omega**2)
	# s2 = -alpha - sqrt(alpha**2 - omega**2)

	# a = np.array([[1,1],[s1, s2]],dtype='float')
	# b = np.array([V0,-((V0+R*Is)/R*C)],dtype='float')
	# x = np.linalg.solve(a,b)
	# A1 = x[1]
	# A2 = x[0]


	A1, A2, s1, s2 = linearSol(alpha,omega,V0, I0, associacao)

	resposta,r,r_degrau = resposta_rlc(alpha, omega,associacao, 0, 0, s1, s2, A1, A2)

	print("Alpha:",alpha)
	print("Omega:",omega)
	print("Raiz s1:",s1)
	print("Raiz s2:",s2)
	print("Amortecimento ",resposta)
	print("Resposta v(t):",r)
	print("Resposta ao degrau v(t):",r_degrau)
	
	tx = np.arange(0,5,0.1)
	f = lambdify(t,r_degrau)
	ty = f(tx)
	plt.plot(tx,ty)
	plt.show()
	print(f(1))


