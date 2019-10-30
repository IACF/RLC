import numpy as np
import math
from sympy import *
import matplotlib.pyplot as plt

m = 10**(-3) ##definicao de mili
degrau = input('Digite (1) para RLC com resposta ao degrau, e (0) para RLC sem resposta ao degrau: ')
R = float(input('Entre com o Resistor: '))
C = float(input('Entre com o Capacitor: '))
L = float(input('Entre com o Indutor: '))
V0 = float(input('Entre com a tensão V(0) do Capacitor: '))
I0 = float(input('Entre com a corrente I(0) do Indutor: '))

associacao = input('Entre com  \"||" para RLC paralelo ou \"--\" para RLC série: ')


Vss = 0
Iss = 0
A1 = symbols('A1')
A2 = symbols('A2')
B1 = symbols('B1')
B2 = symbols('B2')
t = symbols('t')

def linearSol_degrau(alpha, omega, V0, I0, associacao, Vss, Iss): #resolve o sistema linear de acordo com o amortecimento PARA A RESPOSTA AO DEGRAU

	di0 = V0/L # derivada de i(0)
	dv0 = I0/C #derivada de v(0)

	if(omega > alpha):			# subamortecido
		s1 = -alpha
		s2 = sqrt(abs(alpha**2 - omega**2))
		a = np.array([[1,0],[s1, s2]],dtype='float')
		if(associacao == '--'):							# SERIE
			# Tipo de resposta:
			# alpha < omega  --> subamortecido
			# i(t) = e^(-alpha*t)*(A1*cos(wd*t) + A2*sen(wd*t)
			b = np.array([V0-Vss,dv0],dtype='float')
		else:
			b = np.array([I0-Iss,di0],dtype='float')	# PARALELO
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
	elif (omega == alpha):		# critico
		## Para o caso rlc em serie sem fonte
		# s1 = s2 = -alpha = -R/2*L
		s1 = s2 = -alpha
		#s2 = -alpha
		a = np.array([[1,0],[s1, 1]],dtype='float')
		if(associacao == '--'):
			b = np.array([V0-Vss,dv0],dtype='float')
		else:
			b = np.array([I0-Iss,di0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
	else:						# supercritico
		s1 = -alpha + sqrt(alpha**2 - omega**2)
		s2 = -alpha - sqrt(alpha**2 - omega**2)
		a = np.array([[1,1],[s1, s2]],dtype='float')
		if(associacao == '--'):
			b = np.array([V0-Vss,dv0],dtype='float')
		else:
			b = np.array([I0-Iss,di0],dtype='float')
		x = np.linalg.solve(a,b)
		A1 = x[0]
		A2 = x[1]
	return A1, A2, s1, s2

def linearSol(alpha, omega, V0, I0, associacao): #resolve o sistema linear de acordo com o amortecimento

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



def imprime_resultado(r):
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
	
	# plota a resposta
	tx = np.arange(0.,3.5,0.1)
	f = lambdify(t,r) # converte expressão em função
	ty = f(tx)
	plt.xlabel('tempo (s)')
	plt.ylabel(ylabel)
	plt.plot(tx,ty) #plota função
	plt.show()



#----------------------------------------------------------------------------------------#
if associacao == '--':
	while(True):
		try:
			if(degrau == '1'):
				Vss = float(input('Entre com o valor da fonte de Tensão: '))
			break
		except:
			print('Entrada inválida')

	alpha = R/(2*L)
	omega = 1./(sqrt(L*C))
	
	if(degrau == '1'):#verifica se deve chamar a solução pra resposta ao degrau ou pra resposta natural sem fonte
		A1, A2, s1, s2 = linearSol_degrau(alpha,omega,V0, I0, associacao, Vss, 0) # chama função que devolve os coeficientes de acordo com o tipo de amortecimento
	else:
		A1, A2, s1, s2 = linearSol(alpha,omega,V0, I0, associacao) # chama função que devolve os coeficientes de acordo com o tipo de amortecimento
	resposta,r,r_degrau = resposta_rlc(alpha, omega, associacao, s1, s2, A1, A2, Vss, Iss) # devolve a resposta ressonante e a resposta ao degrau do circuito
	
	imprime_resultado(r)

elif associacao == '||':
	
	while(True):
		try:
			if(degrau == '1'):
				Iss = float(input('Entre com o valor da fonte de Corrente : '))
			break
		except:
			print('Entrada inválida')
	
	try: # trata quando a resistência ou capacitância é zero
		alpha = 1/(2*R*C)
	except:
		alpha = 0
	omega = 1/(np.sqrt(L*C))

	if(degrau == '1'): #verifica se deve chamar a solução pra resposta ao degrau ou pra resposta natural sem fonte
		A1, A2, s1, s2 = linearSol_degrau(alpha,omega,V0, I0, associacao, 0, Iss) # chama função que devolve os coeficientes de acordo com o tipo de amortecimento
	else:
		A1, A2, s1, s2 = linearSol(alpha,omega, V0, I0, associacao) # chama função que devolve os coeficientes de acordo com o tipo de amortecimento
		
	resposta,r,r_degrau = resposta_rlc(alpha, omega,associacao, s1, s2, A1, A2, Vss, Iss)

	imprime_resultado(r)


print('\n')
input('Digite enter para finalizar ...')