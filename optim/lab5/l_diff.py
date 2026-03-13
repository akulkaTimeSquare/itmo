import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

# --- 1. Исходные данные ---
A = np.array([[-5, 1], [0, -8]])
b = np.array([[6], [4]])
C = np.array([[1, 0]])
G = np.eye(2)
W = np.array([[8, 3], [3, 2]])
V = 3.0
x0 = np.array([[1], [0]])
x_hat0 = np.array([[0], [0]])

# Оптимальный L (из предыдущего шага)
Q = G @ W @ G.T
P = linalg.solve_continuous_are(A.T, C.T, Q, V)
L_opt = P @ C.T / V
J_opt_theoretical = np.trace(P)

# --- 2. Изменение L (Пункт 3 задания) ---
# Уменьшаем коэффициенты в 2 раза
L_new = L_opt + 1.25

# Проверка устойчивости (собственные числа матрицы A - LC)
A_obs = A - L_new @ C
eig_vals = linalg.eigvals(A_obs)

print("Оптимальная матрица L:\n", L_opt)
print("\nИзмененная матрица L (L_new):\n", L_new)
print("\nСобственные числа наблюдателя с L_new:", eig_vals)
if all(np.real(eig_vals) < 0):
    print("-> Наблюдатель УСТОЙЧИВ.")
else:
    print("-> Наблюдатель НЕУСТОЙЧИВ.")

# --- 3. Моделирование ---
dt = 0.0001
T_end = 15.0
t_values = np.arange(0, T_end, dt)

x = x0.astype(float)
x_hat_opt = x_hat0.astype(float) # Для сравнения (оптимальный)
x_hat_new = x_hat0.astype(float) # Новый (субоптимальный)

# Списки для хранения результатов нового наблюдателя
e_h1_new_list = []
e_h2_new_list = []
J_new_list = []

# Списки для расчета среднего J
J_sum_opt = 0
J_sum_new = 0
count = 0

np.random.seed(42) # Тот же seed для честного сравнения
L_W = np.linalg.cholesky(W * dt)
sqrt_V_dt = np.sqrt(V/dt)

for t in t_values:
    u = np.sin(t)
    
    # Генерация шумов
    noise_proc = G @ (L_W @ np.random.randn(2, 1))
    noise_meas = sqrt_V_dt * np.random.randn()
    
    # Измерение
    y = C @ x + noise_meas
    
    # --- Расчет ошибок для НОВОГО наблюдателя ---
    e_new = x - x_hat_new
    e_h1_new_list.append(e_new[0, 0])
    e_h2_new_list.append(e_new[1, 0])
    J_val_new = float(e_new.T @ e_new)
    J_new_list.append(J_val_new)
    
    # --- Сбор статистики (для сравнения средних) ---
    if t > 2.0: # Считаем среднее после переходного процесса
        e_opt = x - x_hat_opt
        J_sum_opt += float(e_opt.T @ e_opt)
        J_sum_new += J_val_new
        count += 1
    
    # Динамика системы
    dx = A @ x + b * u
    x = x + dx * dt + noise_proc
    
    # Динамика Оптимального наблюдателя
    d_x_hat_opt = A @ x_hat_opt + b * u + L_opt @ (y - C @ x_hat_opt)
    x_hat_opt = x_hat_opt + d_x_hat_opt * dt
    
    # Динамика НОВОГО наблюдателя
    d_x_hat_new = A @ x_hat_new + b * u + L_new @ (y - C @ x_hat_new)
    x_hat_new = x_hat_new + d_x_hat_new * dt

# Расчет средних значений
J_mean_opt = J_sum_opt / count
J_mean_new = J_sum_new / count

print(f"\nСреднее значение J (Оптимальный, симуляция): {J_mean_opt:.4f}")
print(f"Среднее значение J (Измененный L, симуляция): {J_mean_new:.4f}")

# --- 4. Построение графиков ---
plt.figure()
plt.plot(t_values, e_h1_new_list, label='$e_{H1}$', color='green')
plt.title(f'Ошибка $e_{{H1}}$ с измененной матрицей L')
plt.grid(True)
plt.xlabel("t")
plt.ylabel("$e_{H1}(t)$")
plt.legend()
plt.savefig('images/eh1_suboptimal.png')

plt.figure()
plt.plot(t_values, e_h2_new_list, label='$e_{H2}$', color='green')
plt.title(f'Ошибка $e_{{H2}}$ с измененной матрицей L')
plt.xlabel("t")
plt.ylabel("$e_{H2}(t)$")
plt.grid(True)
plt.legend()
plt.savefig('images/eh2_suboptimal.png')

plt.figure()
plt.plot(t_values, J_new_list, label='$J(t)$', color='green')
plt.title('Критерий $J(t)$ с измененной матрицей L')
plt.grid(True)
plt.xlabel("t")
plt.ylabel("$J(t)$")
plt.legend()
plt.savefig('images/j_suboptimal.png')

plt.show()