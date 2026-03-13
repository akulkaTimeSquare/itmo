import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_continuous_are
from scipy.integrate import solve_ivp


def system_dynamics(t, z):
    x = z[:2].reshape(-1, 1)
    
    u = -K @ x
    u_val = u[0, 0]
    
    dx = A @ x + b * u_val
    
    dJ = x.T @ Q @ x + r * (u_val**2)
    
    return [dx[0, 0], dx[1, 0], dJ[0, 0]]

A = np.array([[0, 1],
              [5, -4]])
b = np.array([[1],
              [1]])

x0 = np.array([1, 0])

Q = np.array([[2, 0],
              [0, 1]])
r = 5
R = np.array([[r]])

K = np.array([1.8715, 0.3934]).reshape(1, -1)

t_span = [0, 5]
t_eval = np.linspace(0, 5, 500)
z0 = [x0[0], x0[1], 0]

sol = solve_ivp(system_dynamics, t_span, z0, t_eval=t_eval)

t = sol.t
x1 = sol.y[0]
x2 = sol.y[1]
J_vals = sol.y[2]

u_vals = -K[0,0]*x1 - K[0,1]*x2

P = solve_continuous_are(A, b, Q, R)
J_theoretical = x0.T @ P @ x0

# Графики
plt.figure()
plt.plot(t, x1, label='$x_1(t)$', linewidth=2)
plt.plot(t, x2, label='$x_2(t)$', linewidth=2)
plt.ylabel('$x(t)$')
plt.xlabel('t')
plt.title('Переходные процессы переменных состояния')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("images/x.png")

# График управления
plt.figure()
plt.plot(t, u_vals, label='$u(t)$', linewidth=2)
plt.ylabel('$u(t)$')
plt.xlabel('t')
plt.title('Оптимальное управляющее воздействие')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("images/u.png")

# График критерия
plt.figure()
plt.plot(t, J_vals, label='$J(t)$', linewidth=2)
plt.axhline(y=J_theoretical, color='black', linestyle='--', label=f'Оптимальное $J \\approx {J_theoretical:.2f}$')
plt.ylabel('$J(t)$')
plt.xlabel('t')
plt.title('Накопленное значение критерия качества')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("images/j.png")

plt.show()
