import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

# --- 1. Исходные параметры (базовые) ---
A = np.array([[-5, 1], [0, -8]])
b = np.array([[6], [4]])
C = np.array([[1, 0]])
G = np.eye(2)
V = 3.0
# "Проектное" W (для расчета L)
W_design = np.array([[8, 3], [3, 2]])

# Расчет L_design (как в п.1)
Q_design = G @ W_design @ G.T
P_design = linalg.solve_continuous_are(A.T, C.T, Q_design, V)
L_design = P_design @ C.T / V
J_theoretical_design = np.trace(P_design)

# --- 2. Новые условия (Пункт 4) ---
# Реальное W отличается от проектного (увеличиваем шум)
W_real = W_design * 2.0 

print("Проектное W:\n", W_design)
print("Реальное W (в симуляции):\n", W_real)
print("Матрица L (не меняется):\n", L_design)

# --- 3. Моделирование ---
dt = 0.001
T_end = 10.0
t_values = np.arange(0, T_end, dt)

x = np.array([[1], [0]], dtype=float)
x_hat = np.array([[0], [0]], dtype=float) # Наблюдатель с L_design

e_h1_list = []
e_h2_list = []
J_list = []

J_sum = 0
count = 0

np.random.seed(42)
# ВАЖНО: Генерируем шум на основе W_real
L_W_real = np.linalg.cholesky(W_real * dt)
sqrt_V_dt = np.sqrt(V/dt)

for t in t_values:
    u = np.sin(t)
    
    # Генерация шумов (W_real)
    noise_proc = G @ (L_W_real @ np.random.randn(2, 1))
    noise_meas = sqrt_V_dt * np.random.randn()
    
    y = C @ x + noise_meas
    
    e = x - x_hat
    e_h1_list.append(e[0, 0])
    e_h2_list.append(e[1, 0])
    J_val = float(e.T @ e)
    J_list.append(J_val)
    
    # Сбор статистики
    if t > 2.0:
        J_sum += J_val
        count += 1
    
    # Система: обновляется с W_real
    dx = A @ x + b * u
    x = x + dx * dt + noise_proc
    
    # Наблюдатель: использует L_design
    d_x_hat = A @ x_hat + b * u + L_design @ (y - C @ x_hat)
    x_hat = x_hat + d_x_hat * dt

J_mean_real = J_sum / count
print(f"\nТеоретический J (для проектных условий): {J_theoretical_design:.4f}")
print(f"Средний J в симуляции (с увеличенным шумом W): {J_mean_real:.4f}")

# --- 4. Графики ---
plt.figure()
plt.plot(t_values, e_h1_list, label='$e_{H1}$', color='darkgreen')
plt.title(f'Ошибка $e_{{H1}}$ при увеличенном шуме W')
plt.grid(True)
plt.xlabel("t")
plt.ylabel("$e_{H1}(t)$")
plt.legend()
plt.savefig('images/eh1_w.png')

plt.figure()
plt.plot(t_values, e_h2_list, label='$e_{H2}$', color='darkgreen')
plt.title(f'Ошибка $e_{{H2}}$ при увеличенном шуме W')
plt.grid(True)
plt.xlabel("t")
plt.ylabel("$e_{H2}(t)$")
plt.legend()
plt.savefig('images/eh2_w.png')

plt.figure()
plt.plot(t_values, J_list, label='$J(t)$', color='darkgreen')
plt.title('Критерий $J(t)$ при увеличенном шуме W')
plt.grid(True)
plt.xlabel("t")
plt.ylabel("$J(t)$")
plt.legend()
plt.savefig('images/j_w.png')

plt.show()
