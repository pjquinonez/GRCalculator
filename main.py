import GRCalculator as gr
import sympy as sym
from datetime import datetime

t, r, theta, phi, rs = sym.symbols("t r theta phi rs")
g2_matrix = sym.Matrix([[-(1-(rs/r)),0,0,0],[0,1/(1-(rs/r)),0,0],[0,0,r**2,0],[0,0,0,(r**2)*(sym.sin(theta)**2)]])

g2 = gr.metric(g2_matrix, [t, r, theta, phi])
chris2 = gr.christoffel(g2)
riemann2 = gr.riemann(chris2)
start = datetime.now()
print(riemann2.solve(1,2,1,2))
# for i in range(4):
#     for ii in range(4):
#         for iii in range(4):
#             print(i,ii,iii,":", chris2.solve(i,ii,iii))
print(datetime.now()-start)



