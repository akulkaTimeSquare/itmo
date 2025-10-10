import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

r = 1/8

# Определим систему уравнений
def system(t, x):
    x1, x2 = x
    dx1_dt = -x1
    dx2_dt = r*x2 - x2**3 + x2**5
    
    return [dx1_dt, dx2_dt]

# Создаем сетку начальных условий
x1_vals = np.arange(-1.5, 1.5, 0.2)
x2_vals = np.arange(-1.5, 1.5, 0.2)

# Генерируем фазовые траектории
fig, ax = plt.subplots(figsize=(8, 6))

for x1 in x1_vals:
    for x2 in x2_vals:
        # Решение системы для заданного начального условия
        sol = solve_ivp(system, [0, 10], [x1, x2], t_eval=np.linspace(0, 10, 500))
        ax.plot(sol.y[0], sol.y[1], 'b', alpha=0.6)


ax.scatter([0], [0], label="5")
if r < 1/4:
    ax.scatter([0, 0], [np.sqrt((1 + np.sqrt(1-4*r))/2), -np.sqrt((1 + np.sqrt(1-4*r))/2)], label="1, 2")
    if r > 0:
        ax.scatter([0, 0], [np.sqrt((1 - np.sqrt(1-4*r))/2), -np.sqrt((1 - np.sqrt(1-4*r))/2)], label="3, 4")

# Настройка графика
ax.set_xlabel(r'$x_1$', fontsize=14)
ax.set_ylabel(r'$x_2$', fontsize=14)
ax.set_title('Фазовый портрет маятника', fontsize=16)
ax.legend()
ax.grid()
plt.show()