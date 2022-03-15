import GRCalculator as gr
import sympy as sym
from datetime import datetime

#
# t, x, y, z = sym.symbols("t x y z")
#
# g_matrix = sym.Matrix([[-t,0,0,0],[0,x,0,0],[0,0,y,0],[0,0,0,z]])


# g = gr.metric(g_matrix, [t, x, y, z])
# print(g.get_matrix().inv())
# print(g.get_inv())
# print(g.get_inv_elm(1,1))
# print(sym.diff(g.get_elm(0,0), t))
# chris = gr.christoffel(g)
# print(chris.solve(0,0,0))

t, r, theta, phi, rs = sym.symbols("t r theta phi rs")
g2_matrix = sym.Matrix([[-(1-(rs/r)),0,0,0],[0,1/(1-(rs/r)),0,0],[0,0,r**2,0],[0,0,0,(r**2)*(sym.sin(theta)**2)]])

g2 = gr.metric(g2_matrix, [t, r, theta, phi])
chris2 = gr.christoffel(g2)
print(chris2.solve(2,2,1))



