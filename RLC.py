#import ahkab				# ahkab nao esta funcionando no windows
import numpy as np
import math
from sympy import *
import matplotlib.pyplot as plt

# R = float(input('Entre com o Resistor: '))
# C = float(input('Entre com o Capacitor: '))
# V0 = float(input('Entre com a tensão V(0) do Capacitor: '))
# I0 = float(input('Entre com a corrente I(0) do Indutor: '))
I0 = 0
V0 = 10
# Vs = 24
# L = float(input('Entre com o Indutor: '))
associacao = input('Entre com  \"\\\\" para RLC paralelo ou \"--\" para RLC série: ')

m = 10**(-3) ##definicao de mili
R = 4
C = 0.25
#print('CAPAC:',C)
L = 1

A1 = symbols('A1')
A2 = symbols('A2')
B1 = symbols('B1')
B2 = symbols('B2')
t = symbols('t')

def linearSol(alpha, omega, V0, I0, associacao): #resolve o sistema linear de acordo com o amortecimento

	di0 = -(1/L)*(R*I0+V0) # derivada de i(0)
	dv0 = -(R*I0 + V0)/(R*C) #derivada de v(0)
	print('di0:', di0)
	if(omega > alpha):          # subamortecido
		s1 = -alpha
		s2 = sqrt(abs(alpha**2 - omega**2))
		a = np.array([[1,0],[s1, s2]],dtype='float')
		if(associacao == '--'):							# SERIE
			# Tipo de resposta:
			# alpha < omega  --> subamortecido
			# i(t) = e^(-alpha*t)*(A1*cos(wd*t) + A2*sen(wd*t)
			b = np.array([I0,di0],dtype='float')
		else:
			b = np.array([V0,dv0],dtype='float')		# PARALELO
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
	elif (omega == alpha):          # critico
		## Para o caso rlc em serie sem fonte
		# s1 = s2 = -alpha = -R/2*L
		s1 = s2 = -alpha
		#s2 = -alpha
		a = np.array([[1,0],[s1, 1]],dtype='float')
		b = np.array([V0,dv0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
	else:                            # supercritico
		s1 = -alpha + sqrt(alpha**2 - omega**2)
		s2 = -alpha - sqrt(alpha**2 - omega**2)
		a = np.array([[1,1],[s1, s2]],dtype='float')
		b = np.array([V0,dv0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
	print('A1:',A1, 'A2:',A2)
	return A1, A2, s1, s2





def resposta_rlc(alpha, omega, associacao, Vs, Is, s1, s2, A1, A2): #funcao para verificar  o tipo de resposta e retornar a resposta natural a partir do tipo de amortecimento
	resposta = ""
	
	if alpha > omega:
		resposta = "Amortecimento supercrítico"
		if(associacao == '\\\\'): #resposta para caso seja paralelo
			r = A1*exp(s1*t) + A2*exp(s2*t)
			r_degrau = Is + A1*exp(s1*t) + A2*exp(s2*t)
		elif (associacao == '--') : # resposta caso esteja em série
			r = A1*exp(s1*t) + A2*exp(s2*t)
			r_degrau = Vs + A1*exp(s1*t) + A2*exp(s2*t)#resposta ao DEGRAU para caso seja serie
	elif alpha == omega:
		resposta = "Amortecimento crítico"
		if(associacao == '\\\\'):#resposta para caso seja paralelo
			r = (A1 + A2*t)*exp(-alpha*t)#resposta ressonante para caso seja paralelo
			r_degrau = Is + (A1 + A2*t)*exp(-alpha*t)#resposta ao DEGRAU para caso seja paralelo
		elif (associacao == '--') :
			r = (A2 + A1*t)*exp(-alpha*t)
			r_degrau = Vs + (A2 + A1*t)*exp(-alpha*t)#resposta ao DEGRAU para caso seja serie
	else:
		resposta = "Subamortecimento"
		omega_d = sqrt(abs(omega**2 - alpha**2))
		if(associacao == '\\\\'):
			r = exp(-alpha*t)*(A1*cos(omega_d*t) + A2*sin(omega_d*t)) #resposta ressonante para caso seja paralelo
			r_degrau = Is + exp(-alpha*t)*(A1*cos(omega_d*t) + A2*sin(omega_d*t))#resposta ao DEGRAU para caso seja paralelo
		elif (associacao == '--') :
			r = exp(-alpha*t)*(A1*cos(omega_d*t) + A2*sin(omega_d*t))       # TESTAR
			r_degrau = exp(-alpha*t)*(A1*cos(omega_d*t) + A2*sin(omega_d*t))#resposta ao DEGRAU para caso seja serie
	return resposta,r,r_degrau


def imprime_resultado():
	print("Alpha:",alpha, " Np/s")
	print("Omega:",omega, " rad/s")
	print("Raiz s1:",s1)
	print("Raiz s2:",s2)
	print("########################################")
	print("Tipo de Resposta ",resposta)
	print("Resposta v(t):",r, " V")
	print("Resposta ao degrau v(t):",r_degrau, " V")
	print("########################################")
	# plota a resposta
	tx = np.arange(0,3,0.1)
	f = lambdify(t,r_degrau) # converte expressão em função
	ty = f(tx)
	plt.plot(tx,ty) #plota função
	plt.show()
	print(f(1))


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
	A1, A2, s1, s2 = linearSol(alpha,omega,V0, I0, associacao) # chama função que devolve os coeficientes de acordo com o tipo de amortecimento


	resposta,r,r_degrau = resposta_rlc(alpha, omega,associacao, Vs, 0, s1, s2, A1, A2) # devolve a resposta ressonante e a resposta ao degrau do circuito

	imprime_resultado()	

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


	A1, A2, s1, s2 = linearSol(alpha,omega,V0, I0, associacao) # chama função que devolve os coeficientes de acordo com o tipo de amortecimento

	resposta,r,r_degrau = resposta_rlc(alpha, omega,associacao, 0, 0, s1, s2, A1, A2)

	imprime_resultado()


