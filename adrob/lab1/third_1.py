import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# === Параметры системы ===
a0, a1, b0 = 0, 0, 9
omega_n = 5 + 1/3
aM0 = (28 + 4/9)
aM1 = (10 + 2/3)

# === Матрицы системы и эталонной модели ===
A = np.array([[0, 1],
              [-a0, -a1]])
B = np.array([[0],
              [b0]])
C = np.array([[1, 0]])

A_M = np.array([[0, 1],
                [-aM0, -aM1]])
B_M = np.array([[0],
                [aM0]])
C_M = np.array([[1, 0]])

# === Константы ===
gamma = 1 # скорость адаптации
k = b0 / aM0

theta1 = (-aM0 + a0) / b0
theta2 = (-aM1 + a1) / b0
theta = np.array([theta1, theta2])
print(theta)

# === Псевдообратные коэффициенты ===
b = B.flatten()
b_M = B_M.flatten()

# === Задание входного сигнала ===
def control(t):
    return 9*np.sin(0.2*t) + 9*np.cos(0.1*t) + 15

# === Формирование матриц Ляпунова ===
# Выберем P как решение уравнения Ляпунова: A_M^T P + P A_M = -Q
from scipy.linalg import solve_continuous_lyapunov
Q = np.eye(2)
P = solve_continuous_lyapunov(A_M.T, -Q)

# === Динамика совмещенной системы (объект + эталонная модель + адаптация) ===
def system(t, z):
    # Вектор состояния z = [x(2), xM(2), theta_hat(2)]
    x = z[0:2]
    xM = z[2:4]
    theta_hat = z[4:6]

    g = control(t)

    # --- Управление ---
    u = theta_hat @ x + (1/k) * g

    # --- Объект и эталонная модель ---
    dx = (A @ x + B.flatten() * u)
    dxM = (A_M @ xM + B_M.flatten() * g)

    # --- Ошибка и адаптация ---
    e = xM - x
    dtheta = gamma * np.outer(x, b.T @ P @ e)  # γ x b^T P e
    dtheta = dtheta.flatten()[:2]

    return np.hstack((dx, dxM, dtheta))

# === Начальные условия ===
x0 = np.array([1.0, 1.0])
xM0 = np.array([0.0, 0.0])
theta0 = np.ones(2)
z0 = np.hstack((x0, xM0, theta0))

# === Интегрирование ===
t_span = [0, 350]
dt = 0.01
t_eval = np.arange(t_span[0], t_span[1], dt)
sol = solve_ivp(system, t_span, z0, t_eval=t_eval, rtol=1e-6, atol=1e-8)

x = sol.y[0:2]
xM = sol.y[2:4]
theta_hat = sol.y[4:6]
e = xM - x

# === Визуализация ===
plt.figure(figsize=(8, 5))
plt.plot(sol.t, x[0], label='$x_1(t)$', color='blue', linewidth=2)
plt.plot(sol.t, x[1], label='$x_2(t)$', color='green', linewidth=2)
plt.plot(sol.t, xM[0], '--', label='$x_{M1}(t)$', color='red', linewidth=2)
plt.plot(sol.t, xM[1], '--', label='$x_{M2}(t)$', color='black', linewidth=2)
plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('Состояния объекта и эталонной модели при расчетных параметрах')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.xlim(t_span[0], t_span[1])
plt.savefig('images/3_known_params_response.png', dpi=450)
plt.show()

# === Ошибка ===
plt.figure(figsize=(8, 5))
plt.plot(sol.t, e[0], label='$e_1(t)$', color='red', linewidth=2)
plt.plot(sol.t, e[1], label='$e_2(t)$', color='black', linewidth=2)
plt.axhline(0, color='gray', linestyle='--', linewidth=1)
plt.xlabel('t')
plt.ylabel('e(t)')
plt.title('Ошибка между объектом и эталонной моделью при расчетных параметрах')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.xlim(t_span[0], t_span[1])
plt.savefig('images/3_known_params_error.png', dpi=450)
plt.show()

# === Параметры θ̂ ===
plt.figure(figsize=(8, 5))
plt.plot(sol.t,  theta[0]*np.ones(np.size(theta_hat[0])) - theta_hat[0], label=r'$\tilde{\theta}_1(t)$', linewidth=2)
plt.plot(sol.t, theta[1]*np.ones(np.size(theta_hat[1])) - theta_hat[1], label=r'$\tilde{\theta}_2(t)$', linewidth=2)
plt.xlabel('t')
plt.ylabel(r'$\tilde{\theta}(t)$')
plt.title(r'Ошибка оценки $\tilde{\theta}(t)$ при расчетных параметрах')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.xlim(t_span[0], t_span[1])
plt.savefig('images/3_known_params_estimation.png', dpi=450)
plt.show()