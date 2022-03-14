import GRCalculator as gr
import sympy as sym


t, x, y, z = sym.symbols("t x y z")

g_matrix = sym.Matrix([[-t,0,0,0],[0,x,0,0],[0,0,y,0],[0,0,0,z]])


g = gr.metric(g_matrix, [t, x, y, z])
# print(g.get_matrix().inv())
# print(g.get_inverse())
# print(g.get_inverse_element(1,1))
# print(sym.diff(g.get_element(0,0), t))
chris = gr.christoffel(g)
print(chris.solve(0,0,0))

