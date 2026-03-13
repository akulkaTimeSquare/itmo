import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

# --- System Matrices ---
A = np.array([[-5, 1], [0, -8]])
b = np.array([[6], [4]])
C = np.array([[1, 0]])
G = np.eye(2)
W = np.array([[8, 3], [3, 2]])
V = 3.0

# --- LQR Parameters (from prompt) ---
Q_reg = np.array([[2, 0], [0, 1]])
r_reg = 5.0

# --- 1. Calculate Regulator Gain K ---
# Solve ARE for Control: A^T P + P A - P b r^-1 b^T P + Q = 0
# scipy solve_continuous_are solves: A^T X + X A - X B R^-1 B^T X + Q = 0
# Mapping for Regulator:
# X -> P_reg
# A -> A
# B -> b
# R -> r_reg
# Q -> Q_reg
P_reg = linalg.solve_continuous_are(A, b, Q_reg, r_reg)
K_reg = (1/r_reg) * b.T @ P_reg

print("Regulator P:\n", P_reg)
print("Regulator K:\n", K_reg)

# --- 2. Calculate Observer Gain L (Standard) ---
Q_obs = G @ W @ G.T
P_obs = linalg.solve_continuous_are(A.T, C.T, Q_obs, V)
L_obs = P_obs @ C.T / V

print("Observer L:\n", L_obs)

# --- 3. LQG Simulation ---
dt = 0.001
T_end = 20.0
t_values = np.arange(0, T_end, dt)

x = np.array([[10], [-10]], dtype=float) # x(0)
x_hat = np.array([[0], [0]], dtype=float) # x_hat(0)

x_hist = []
x_hat_hist = []
u_hist = []
J_lqg_list = [] # Just to see

np.random.seed(42)
L_W = np.linalg.cholesky(W * dt)
sqrt_V_dt = np.sqrt(V/dt)

for t in t_values:
    # 1. Calculate Control u = - K * x_hat
    u_val = - K_reg @ x_hat
    u = float(u_val) # scalar
    
    # 2. Store history
    x_hist.append(x.flatten())
    x_hat_hist.append(x_hat.flatten())
    u_hist.append(u)
    
    # 3. Noise
    noise_proc = G @ (L_W @ np.random.randn(2, 1))
    noise_meas = sqrt_V_dt * np.random.randn()
    
    # 4. Measurement
    y = C @ x + noise_meas
    
    # 5. Dynamics
    # System: dx = Ax + bu + w
    dx = A @ x + b * u
    x = x + dx * dt + noise_proc
    
    # Observer: d_x_hat = A x_hat + b u + L(y - C x_hat)
    d_x_hat = A @ x_hat + b * u + L_obs @ (y - C @ x_hat)
    x_hat = x_hat + d_x_hat * dt

x_hist = np.array(x_hist)
x_hat_hist = np.array(x_hat_hist)
u_hist = np.array(u_hist)

# --- 4. Plots ---
# State x1
plt.figure()
plt.plot(t_values, x_hist[:, 0], label='$x_1$', linewidth=2)
plt.plot(t_values, x_hat_hist[:, 0], '--', label='Оценка $\hat{x}_1$', linewidth=2)
plt.title('Состояние $x_1$ и оценка $\hat{x}_1$')
plt.grid(True)
plt.xlabel("t")
plt.ylabel("f(t)")
plt.legend()
plt.savefig('images/x1_lqg.png')

# State x2
plt.figure()
plt.plot(t_values, x_hist[:, 1], label='$x_2$', linewidth=2)
plt.plot(t_values, x_hat_hist[:, 1], '--', label='Оценка $\hat{x}_2$', linewidth=2)
plt.title('Состояние $x_2$ и оценка $\hat{x}_2$')
plt.grid(True)
plt.xlabel("t")
plt.ylabel("f(t)")
plt.legend()
plt.savefig('images/x2_lqg.png')

# Control u
plt.figure()
plt.plot(t_values, u_hist, label='u(t)')
plt.title('Управляющий сигнал $u(t) = -K \hat{x}$')
plt.grid(True)
plt.xlabel("t")
plt.ylabel("u(t)")
plt.legend()
plt.savefig('images/u_lqg.png')

plt.show()