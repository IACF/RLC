import ahkab
import matplotlib
import six
import Tkinter

mycir = ahkab.Circuit('Simple Example Circuit')

mycir.add_resistor('R1', 'n1', mycir.gnd, value=5)
mycir.add_vsource('V1', 'n2', 'n1', dc_value=8)
mycir.add_resistor('R2', 'n2', mycir.gnd, value=2)
mycir.add_vsource('V2', 'n3', 'n2', dc_value=4)
mycir.add_resistor('R3', 'n3', mycir.gnd, value=4)
mycir.add_resistor('R4', 'n3', 'n4', value=1)
mycir.add_vsource('V3', 'n4', mycir.gnd, dc_value=10)
mycir.add_resistor('R5', 'n2', 'n4', value=4)

opa = ahkab.new_op()
r = ahkab.run(mycir, opa)['op']

print (r)