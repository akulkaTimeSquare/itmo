import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# === Заданные параметры ===
a0, a1, b0 = 0, 0, 9
omega_n = 5 + 1/3
aM0 = (28 + 4/9)
aM1 = (10 + 2/3)

# === Матрицы системы ===
A = np.array([[0, 1],
              [-a0 + 1, -a1 + 1]])
B = np.array([[0],
              [b0]])
C = np.array([[1, 0]])

A_M = np.array([[0, 1],
                [-aM0, -aM1]])
B_M = np.array([[0],
                [aM0]])
C_M = np.array([[1, 0]])

# === Расчёт параметров управления ===
theta1 = (-aM0 + a0) / b0
theta2 = (-aM1 + a1) / b0
theta = np.array([[theta1, theta2]])
print(theta)
b_M = B_M[1, 0]
k = b0 / aM0
print(k)
b = k * b_M

def control(t):
    return 9*np.sin(0.2*t) + 9*np.cos(0.1*t) + 15

def model_ref(t, xM):
    g = control(t)
    dxM = (A_M @ xM) + (B_M.flatten() * g)
    return dxM

def model_obj(t, x):
    g = control(t)
    u = float(theta @ x + (1 / k) * g)
    dx = (A @ x) + (B.flatten() * u)
    return dx

# === Интегрирование ===
t_span = [0, 10]
t_eval = np.linspace(t_span[0], t_span[1], 1000)
xM0 = np.array([0, 0])
x0 = np.array([1, 1])

sol_M = solve_ivp(model_ref, t_span, xM0, t_eval=t_eval)
sol_O = solve_ivp(model_obj, t_span, x0, t_eval=t_eval)

yM = (C_M @ sol_M.y).flatten()
y = (C @ sol_O.y).flatten()
e = sol_M.y - sol_O.y

# === График отклика ===
plt.figure(figsize=(8, 5))
plt.plot(t_eval, sol_O.y[0], label='$x_1(t)$', color='blue', linewidth=2)
plt.plot(t_eval, sol_O.y[1], label='$x_2(t)$', color='green', linewidth=2)
plt.plot(t_eval, sol_M.y[0], label='$x_{M1}(t)$', color='red', linewidth=2, linestyle='--')
plt.plot(t_eval, sol_M.y[1], label='$x_{M2}(t)$', color='black', linewidth=2, linestyle='--')
plt.xlim(t_span[0], t_span[1])
plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('Состояния объекта и эталонной модели')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('images/unknown_params_response_stable.png', dpi=450)
plt.show()

# === График ошибки ===
plt.figure(figsize=(8, 5))
plt.plot(t_eval, e[0], color='red', linewidth=2, label='$e_1(t)$')
plt.plot(t_eval, e[1], color='black', linewidth=2, label='$e_2(t)$')
plt.axhline(0, color='gray', linestyle='--', linewidth=1.5)
plt.xlim(t_span[0], t_span[1])
plt.xlabel('t')
plt.ylabel('e(t)')
plt.title('Ошибка между эталонной моделью и объектом')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('images/unknown_params_error_stable.png', dpi=450)
plt.show()