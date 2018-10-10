import numpy as np
import math
from sympy import *
import matplotlib.pyplot as plt

Vss = 0
Iss = 0
A1 = symbols('A1')
A2 = symbols('A2')
B1 = symbols('B1')
B2 = symbols('B2')
t = symbols('t')

def linearSol(alpha, omega, V0, I0, associacao, R, C, L): #resolve o sistema linear de acordo com o amortecimento

	di0 = -(1/L)*(R*I0+V0) # derivada de i(0)
	try:			# Trata para o caso da resistencia ou capacitancia ser zero
		dv0 = -(R*I0 + V0)/(R*C) #derivada de v(0)
	except:
		dv0 = 0

	if(omega > alpha):				# subamortecido
		s1 = -alpha
		s2 = sqrt(abs(alpha**2 - omega**2))
		a = np.array([[1,0],[s1, s2]],dtype='float')
		if(associacao == '--'):
			# Tipo de resposta:
			# alpha < omega  --> subamortecido
			# i(t) = e^(-alpha*t)*(A1*cos(wd*t) + A2*sen(wd*t)
			b = np.array([I0,di0],dtype='float')
		else:
			b = np.array([V0,dv0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
	elif (omega == alpha):			 # critico
		## Para o caso rlc em serie sem fonte
		s1 = s2 = -alpha
		a = np.array([[1,0],[s1, 1]],dtype='float')
		b = np.array([V0,dv0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]		
	else:							# supercritico
		s1 = -alpha + sqrt(alpha**2 - omega**2)
		s2 = -alpha - sqrt(alpha**2 - omega**2)
		a = np.array([[1,1],[s1, s2]],dtype='float')
		b = np.array([V0,dv0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
	return A1, A2, s1, s2



def resposta_rlc(alpha, omega, associacao, s1, s2, A1, A2, Vss, Iss): #funcao para verificar  o tipo de resposta e retornar a resposta natural a partir do tipo de amortecimento
	resposta = ""
	print("S1: ", s1, "S2:", s2)
	if alpha > omega:
		resposta = "Amortecimento supercrítico"
		if(associacao == '||'): #resposta para caso seja paralelo
			r = A1*exp(s1*t) + A2*exp(s2*t)
			r_degrau = Iss + A1*exp(s1*t) + A2*exp(s2*t)
		elif (associacao == '--') : # resposta caso esteja em série
			r = A1*exp(s1*t) + A2*exp(s2*t)
			r_degrau = Vss + A1*exp(s1*t) + A2*exp(s2*t)#resposta ao DEGRAU para caso seja serie
	elif alpha == omega:
		resposta = "Amortecimento crítico"
		if(associacao == '||'):#resposta para caso seja paralelo
			r = (A1 + A2*t)*exp(-alpha*t)#resposta ressonante para caso seja paralelo
			r_degrau = Iss + (A1 + A2*t)*exp(-alpha*t)#resposta ao DEGRAU para caso seja paralelo
		elif (associacao == '--') :
			r = (A2 + A1*t)*exp(-alpha*t)
			r_degrau = Vss + (A1 + A2*t)*exp(-alpha*t)#resposta ao DEGRAU para caso seja serie
	else:
		resposta = "Subamortecimento"
		omega_d = sqrt(abs(omega**2 - alpha**2))
		if(associacao == '||'):
			r = exp(-alpha*t)*(A1*cos(omega_d*t) + A2*sin(omega_d*t)) #resposta ressonante para caso seja paralelo
			r_degrau = Iss + exp(-alpha*t)*(A1*cos(omega_d*t) + A2*sin(omega_d*t))#resposta ao DEGRAU para caso seja paralelo
		elif (associacao == '--') :
			r = exp(-alpha*t)*(A1*cos(omega_d*t) + A2*sin(omega_d*t))
			r_degrau = Vss + exp(-alpha*t)*(A1*cos(omega_d*t) + A2*sin(omega_d*t))#resposta ao DEGRAU para caso seja serie
	return resposta,r,r_degrau



def imprime_resultado(r, resposta, alpha, omega, I0, V0, L, degrau, r_degrau, associacao):
	print("Alpha:",alpha, " Np/s")
	print("Omega:",omega, " rad/s")
	print("########################################")
	print("Tipo de Resposta ",resposta)
	if(I0 == 0):
		print("Resposta i(t):",r, " A")
		ylabel = 'corrente (A)'				# Titulo da Ordenada do grafico
		if(degrau == '1'):
			print("Resposta ao degrau i(t):",r_degrau, " A")
		print("########################################")
	elif(V0 == 0):
		if(alpha == 0):   # restriçao para um caso em especifico
			r = r.diff(t)*(-L)
		print("Resposta v(t):",r, " V")
		ylabel = 'tensão (V)'
		if(degrau == '1'):
			print("Resposta ao degrau v(t):",r_degrau, " V")
		print("########################################")
	elif(associacao == '--'):
		print("Resposta i(t):", r, " A")
		ylabel = 'corrente (A)'				# Titulo da Ordenada do grafico
		if(degrau == '1'):
			print("Resposta ao degrau i(t):", r_degrau, " A")
		print("########################################")
	elif(associacao == '||'):
		print("Resposta v(t):",r, " V")
		ylabel = 'tensão (V)'
		if(degrau == '1'):
			print("Resposta ao degrau v(t):",r_degrau, " V")
		print("########################################")
	
	# # plota a resposta
	# tx = np.arange(0.,3.5,0.1)
	# f = lambdify(t,r) # converte expressão em função
	# ty = f(tx)
	# plt.xlabel('tempo (s)')
	# plt.ylabel(ylabel)
	# plt.plot(tx,ty) #plota função
	# plt.show()