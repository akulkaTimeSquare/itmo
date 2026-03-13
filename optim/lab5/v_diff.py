import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

# --- 1. Исходные параметры (как в п.1) ---
A = np.array([[-5, 1], [0, -8]])
b = np.array([[6], [4]])
C = np.array([[1, 0]])
G = np.eye(2)
W = np.array([[8, 3], [3, 2]])
V_design = 3.0 # Проектное значение

# Рассчитываем L для проектных условий
Q = G @ W @ G.T
P_design = linalg.solve_continuous_are(A.T, C.T, Q, V_design)
L_design = P_design @ C.T / V_design
J_theoretical = np.trace(P_design)

# --- 2. Новые условия (Пункт 5) ---
# Реальное V значительно больше проектного
V_real = 12.0 

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
L_W = np.linalg.cholesky(W * dt)
# Шум измерений генерируем с V_real
sqrt_V_real_dt = np.sqrt(V_real/dt)

for t in t_values:
    u = np.sin(t)
    
    noise_proc = G @ (L_W @ np.random.randn(2, 1))
    noise_meas = sqrt_V_real_dt * np.random.randn()
    
    y = C @ x + noise_meas
    
    e = x - x_hat
    e_h1_list.append(e[0, 0])
    e_h2_list.append(e[1, 0])
    J_val = float(e.T @ e)
    J_list.append(J_val)
    
    if t > 2.0:
        J_sum += J_val
        count += 1
    
    dx = A @ x + b * u
    x = x + dx * dt + noise_proc
    
    # Наблюдатель использует "старый" L_design
    d_x_hat = A @ x_hat + b * u + L_design @ (y - C @ x_hat)
    x_hat = x_hat + d_x_hat * dt

J_mean_real = J_sum / count

print(f"L (проектный): \n{L_design}")
print(f"Теоретический J (при V=3): {J_theoretical:.4f}")
print(f"Реальный J (при V=12): {J_mean_real:.4f}")

# --- 4. Графики ---
plt.figure()
plt.plot(t_values, e_h1_list, label='$e_{H1}$', color='darkred')
plt.title(f'Ошибка $e_{{H1}}$ при росте шума измерений')
plt.grid(True)
plt.xlabel("t")
plt.ylabel("$e_{H1}(t)$")
plt.legend()
plt.savefig('images/eh1_v.png')

plt.figure()
plt.plot(t_values, e_h2_list, label='$e_{H2}$', color='darkred')
plt.title(f'Ошибка $e_{{H2}}$ при росте шума измерений')
plt.grid(True)
plt.xlabel("t")
plt.ylabel("$e_{H2}(t)$")
plt.legend()
plt.savefig('images/eh2_v.png')

plt.figure()
plt.plot(t_values, J_list, label='$J(t)$', color='darkred')
plt.title('Критерий $J(t)$ при росте шума измерений')
plt.grid(True)
plt.xlabel("t")
plt.ylabel("$J(t)$")
plt.legend()
plt.savefig('images/j_v.png')

plt.show()
