import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

# Matrices
A = np.array([[-5, 1], [0, -8]])
b = np.array([[6], [4]])
C = np.array([[1, 0]])
G = np.eye(2)
W = np.array([[8, 3], [3, 2]])
V = 3.0

# Initial conditions
x0 = np.array([[1], [0]])
x_hat0 = np.array([[0], [0]])

Q_scipy = G @ W @ G.T

P = linalg.solve_continuous_are(A.T, C.T, Q_scipy, V)
print("Matrix P:")
print(P)

V_inv = 1/V
L = P @ C.T * V_inv
print("\nMatrix L:")
print(L)

# Theoretical J
J_theoretical = np.trace(P)
print(f"\nTheoretical J (Trace of P): {J_theoretical}")

# 2. Simulation
dt = 0.0001
T_end = 15.0
time_steps = int(T_end / dt)
t_values = np.linspace(0, T_end, time_steps)

x = x0.astype(float)
x_hat = x_hat0.astype(float)

e_h1_list = []
e_h2_list = []
J_list = [] # J(t) = e^T e

np.random.seed(42)

for t in t_values:
    # Input
    u = np.sin(t)
    
    L_W = np.linalg.cholesky(W * dt)
    noise_proc = G @ (L_W @ np.random.randn(2, 1))
    
    # Measurement noise
    # v ~ N(0, V/dt)
    # std = sqrt(V/dt)
    noise_meas = np.sqrt(V/dt) * np.random.randn()
    
    # Measurement
    y = C @ x + noise_meas
    
    # Store errors
    e = x - x_hat
    e_h1_list.append(e[0, 0])
    e_h2_list.append(e[1, 0])
    J_val = float(e.T @ e)
    J_list.append(J_val)
    
    # Update system
    dx = A @ x + b * u
    x = x + dx * dt + noise_proc
    
    # Update observer
    # d_x_hat = A x_hat + b u + L (y - C x_hat)
    d_x_hat = A @ x_hat + b * u + L @ (y - C @ x_hat)
    x_hat = x_hat + d_x_hat * dt

# Convert to arrays
e_h1 = np.array(e_h1_list)
e_h2 = np.array(e_h2_list)
J_arr = np.array(J_list)
print(np.mean(J_arr))

# Plotting
plt.figure()
plt.plot(t_values, e_h1, label='$e_{H1}$')
plt.title('Ошибка $e_{H1}(t)$')
plt.grid(True)
plt.xlabel("t")
plt.ylabel("$e_{H1}(t)$")
plt.legend()
plt.savefig('images/eh1.png')

plt.figure()
plt.plot(t_values, e_h2, label='$e_{H2}$')
plt.title('Ошибка $e_{H2}(t)$')
plt.grid(True)
plt.xlabel("t")
plt.ylabel("$e_{H2}(t)$")
plt.legend()
plt.savefig('images/eh2.png')

plt.figure()
plt.plot(t_values, J_arr, label='$J(t) = e_H^T(t) e_H(t)$')
plt.title('Критерий качества $J(t)$')
plt.grid(True)
plt.xlabel("t")
plt.ylabel("$J(t)$")
plt.legend()
plt.savefig('images/j.png')

plt.show()
