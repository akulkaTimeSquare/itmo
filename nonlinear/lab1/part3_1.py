# Python simulation and LQR design (visible to user)
import numpy as np
from scipy.linalg import solve_continuous_are, eig
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd
import math

# Given linearization
A = np.array([[-1, 1],
              [-1, -1]])
B = np.array([[1, 0],
              [0, 3]])

# LQR design (continuous-time)
Q = np.diag([1.0, 1.0])   # state penalty (tuneable)
R = np.diag([1.0, 1.0])     # input penalty (tuneable)

# Solve continuous-time Algebraic Riccati Equation
P = solve_continuous_are(A, B, Q, R)
# Compute gain K = R^-1 B^T P
K = np.linalg.inv(R) @ B.T @ P

# Closed-loop linear eigenvalues
Acl = A - B @ K
eigvals, _ = eig(Acl)

print("LQR gain K:")
print(np.round(K, 4))
print("\nClosed-loop eigenvalues (linearized):")
print(np.round(eigvals, 4))

# Nonlinear system dynamics
def f_nl(t, x):
    # state-feedback u = -K x
    u = -K @ x
    u1, u2 = u[0], u[1]
    dx1 = -x[0] + 2*(x[0]**3) + x[1] + math.sin(u1)
    dx2 = -x[0] - x[1] + 3*math.sin(u2)
    return [dx1, dx2]

# Simulation settings
t_span = (0.0, 10.0)
t_eval = np.linspace(t_span[0], t_span[1], 1001)

# Several initial conditions near origin (local stabilization)
inits = [
    [0.5, 0.4],
    [-0.2, -0.3]
]

results = []
for x0 in inits:
    sol = solve_ivp(f_nl, t_span, x0, t_eval=t_eval, rtol=1e-7, atol=1e-9)
    results.append((x0, sol))

# Plot states vs time for each initial condition
for (x0, sol) in results:
    plt.figure()
    plt.plot(sol.t, sol.y[0], label='x1')
    plt.plot(sol.t, sol.y[1], label='x2')
    plt.xlabel("t")
    plt.ylabel("x(t)")
    plt.legend(loc='best')
    plt.grid(True)
    plt.tight_layout(pad=1, w_pad=1, h_pad=1)
    plt.savefig(f'images/part3_1_{x0}.png', dpi=600)
    plt.show()

# Phase plot (x1 vs x2) for all trajectories
plt.figure()
for (x0, sol) in results:
    plt.plot(sol.y[0], sol.y[1], label=f"x0={x0}")
    plt.scatter([x0[0]], [x0[1]])  # initial point
plt.xlabel("x1")
plt.ylabel("x2")
plt.legend(loc="best")
plt.grid(True)
plt.plot(0,0,'ro')
plt.tight_layout(pad=1, w_pad=1, h_pad=1)
plt.savefig(f'images/part3_1_phase.png', dpi=600)
plt.show()

# Also plot control inputs u(t) for each run
for (x0, sol) in results:
    U = -K @ sol.y
    plt.figure()
    plt.plot(sol.t, U[0,:], label='u1')
    plt.plot(sol.t, U[1,:], label='u2')
    plt.xlabel("t")
    plt.ylabel("u(t)")
    plt.legend(loc='best')
    plt.grid(True)
    plt.tight_layout(pad=1, w_pad=1, h_pad=1)
    plt.savefig(f'images/part3_1_u_{x0}.png', dpi=600)
    plt.show()