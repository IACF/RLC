import ahkab
import numpy as np

R = float(input('Entre com o Resistor:\n'))
C = float(input('Entre com o Capacitor:\n'))
L = float(input('Entre com o Indutor:\n'))

associacao = input('Entre com  \"\\\\" para RLC paralelo ou \"--\" para RLC série: ')
if associacao == '--':
	
	bpf = ahkab.Circuit('RLC Série sem fonte')
	bpf.add_inductor('L1', 'in', 'n1', L)
	bpf.add_capacitor('C1', 'n1', 'out', C)
	bpf.add_resistor('R1', 'out', bpf.gnd, R)
	# we also give V1 an AC value since we wish to run an AC simulation
	# in the following
	# bpf.add_vsource('V1', 'in', bpf.gnd, dc_value=1, ac_value=1)

	print(bpf)
	omega = 1./(np.sqrt(L*C))
	print ('Frequência Natural não amortecida: %g rad/s' %omega)

	alfa = R/(2*L)

	print ('Frequência Neper: %g Np/s' %alfa)

	f0 = 1./(2*np.pi*np.sqrt(L*C))
	print ('Frequência de ressonância: %g Hz' %f0)

	if alfa > omega:
		print('A Resposta natural é de amortecimento supercrítico ')
	elif alfa == omega:
		print('A Resposta natural é de amortecimento crítico ')
	else:
		print('A Resposta natural é de subamortecimento ')
elif associacao == '\\\\':
	

	bpf = ahkab.Circuit('RLC Paralelo sem fonte')
	bpf.add_inductor('L1', 'n', bpf.gnd, L)
	bpf.add_capacitor('C1', 'n', bpf.gnd, C)
	bpf.add_resistor('R1', 'n', bpf.gnd, R)
	# we also give V1 an AC value since we wish to run an AC simulation
	# in the following
	# bpf.add_vsource('V1', 'in', bpf.gnd, dc_value=1, ac_value=1)

	print(bpf)
	omega = 1./(np.sqrt(L*C))
	print ('Frequência Natural não amortecida: %g rad/s' %omega)

	alfa = 1./(2*R*C)

	print ('Frequência Neper: %g Np/s' %alfa)

	f0 = 1./(2*np.pi*np.sqrt(L*C))
	print ('Frequência de ressonância: %g Hz' %f0)

	if alfa > omega:
		print('A Resposta natural é de amortecimento supercrítico ')
	elif alfa == omega:
		print('A Resposta natural é de amortecimento crítico ')
	else:
		print('A Resposta natural é de subamortecimento ')



