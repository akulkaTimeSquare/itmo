import numpy as np
from scipy import linalg
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def system_dynamics(t, state):
    x1, x2 = state[0], state[1]
    xm1, xm2 = state[2], state[3]
    theta_hat = state[4:6]
    
    g = 28 * np.sin(t) + 26
    
    R = np.array([np.sin(x1), np.cos(x1)])
    
    am_T_x = 8*x1 + 6*x2
    
    u = 8*g - am_T_x - np.dot(theta_hat, R)
    
    dx1 = x2
    dx2 = np.dot(theta_true, R) + u
    
    dxm1 = xm2
    dxm2 = -8*xm1 - 6*xm2 + 8*g
    
    e = np.array([xm1 - x1, xm2 - x2])
    
    Pe = P @ e
    e_P_E2 = Pe[1]
    d_theta_hat = - Gamma @ R * e_P_E2
    
    return [dx1, dx2, dxm1, dxm2, d_theta_hat[0], d_theta_hat[1]]


theta_true = np.array([-10, 2.0])
theta_hat_0 = np.array([0.5, 0.5])

gamma1 = gamma2 = 100
Gamma = np.diag([gamma1, gamma2])

Am = np.array([[0, 1],
               [-8, -6]])
Bm = np.array([0, 8])

Q = np.eye(2)
P = linalg.solve_continuous_lyapunov(Am.T, -Q)

# x(0)
x0 = [0, 0]
xm0 = [0, 0]
initial_state = x0 + xm0 + list(theta_hat_0)

# Интегрирование
t_span = (0, 75)
t_eval = np.arange(0, 75, 0.001)
sol = solve_ivp(system_dynamics, t_span, initial_state, t_eval=t_eval)

t = sol.t
e1 = sol.y[2] - sol.y[0]
e2 = sol.y[3] - sol.y[1]
theta_hat1 = sol.y[4]
theta_hat2 = sol.y[5]

u_vals = []
for i in range(len(t)):
    x1_val = sol.y[0][i]
    x2_val = sol.y[1][i]
    g_val = 28 * np.sin(t[i]) + 26
    w_val = np.array([np.sin(x1_val), np.cos(x1_val)])
    theta_h_val = np.array([theta_hat1[i], theta_hat2[i]])
    u = 8*g_val - (8*x1_val + 6*x2_val) - np.dot(theta_h_val, w_val)
    u_vals.append(u)

# Графики
fig, axs = plt.subplots(3, 1, figsize=(12, 9))

axs[0].plot(t, e1, label='e1 = xm1 - x1')
axs[0].plot(t, e2, label='e2 = xm2 - x2')
axs[0].set_title('Ошибка слежения e(t)')
axs[0].grid(True)
axs[0].legend()

axs[1].plot(t, u_vals, label='Управление u(t)', color='orange')
axs[1].set_title('Сигнал управления u(t)')
axs[1].grid(True)
axs[1].legend()

axs[2].plot(t, theta_hat1, label='Оценка theta1')
axs[2].plot(t, theta_hat2, label='Оценка theta2')
axs[2].axhline(y=theta_true[0], color='r', linestyle='--', label='Истинная theta1')
axs[2].axhline(y=theta_true[1], color='g', linestyle='--', label='Истинная theta2')
axs[2].set_title('Адаптация параметров')
axs[2].grid(True)
axs[2].legend(loc='lower right')

plt.tight_layout()
plt.show()
