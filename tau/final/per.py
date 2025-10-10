import sympy as sp
import numpy as np

n = 11
rng = np.random.Generator(np.random.Philox(n))
M = rng.integers(100000, 1000000) / 1000 / np.sqrt(2)
m = rng.integers(1000, 10000) / 1000 * np.sqrt(3)
l = rng.integers(100, 1000) / np.sqrt(5) / 100
g = 9.81


A = np.array([[0, 1, 0, 0], 
              [0, 0, 3*m*g / (4*M+m), 0], 
              [0, 0, 0, 1], 
              [0, 0, 6*(M+m)*g/l/(4*M+m), 0]])

B = np.array([[0], 
              [4 / (4*M+m)],
              [0], 
              [6/l/(4*M+m)]])

C = np.array([[1, 0, 0, 0],
              [0, 0, 1, 0]])

D = np.array([[0], 
              [6/l/(4*M+m)], 
              [0], 
              [12*(M+m) / m / l**2 / (4*M+m)]])


s = sp.symbols('s')

I = sp.eye(4)
sIA = s * I - A
sIA_inv = sIA.inv()
W = C @ sIA_inv @ B

sp.init_printing() 
print("Wuy(s) = C(sI - A)^{-1}B:\n")
sp.pprint(W.applyfunc(sp.together).evalf(4))

s = sp.symbols('s')

I = sp.eye(4)
sIA = s * I - A
sIA_inv = sIA.inv()
W = C @ sIA_inv @ D

sp.init_printing() 
print("\nWfy(s) = C(sI - A)^{-1}D:\n")
sp.pprint(W.applyfunc(sp.together).evalf(4))     