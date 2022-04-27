import GRCalculator as gr
import sympy as sym
from datetime import datetime

t, r, theta, phi, m, G = sym.symbols("t r theta phi m G")
g_matrix = sym.Matrix([[-(1-((2*G*m)/r)),0,0,0],[0,1/(1-((2*G*m)/r)),0,0],[0,0,r**2,0],[0,0,0,(r**2)*(sym.sin(theta)**2)]])
g = gr.metric(g_matrix, [t, r, theta, phi])

chris = gr.christoffel(g)
riemann = gr.riemann(chris)

start = datetime.now()
print(sym.simplify(chris.solve(0,1,0)))
print(sym.simplify(riemann.solve(1,2,1,2)))
print(datetime.now()-start)


