import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# Требуемые условия
zeta = 1.0
t_n = 0.9
y_target = 0.95

# Определим уравнение для критического затухания
def func(wn):
    return 1 - (1 + wn * t_n) * np.exp(-wn * t_n) - y_target

# Решаем уравнение для wn
wn_solution = 4.8/t_n  # Начальное приближение

# Параметры модели
aM0 = wn_solution**2
aM1 = 2*zeta*wn_solution

# Строим график переходной функции
t = np.linspace(0, 2.5, 1001)
y = 1 - (1 + wn_solution * t) * np.exp(-wn_solution * t)

# Границы ±5%
y_upper, y_lower = 1.05, 0.95

# Построение графика
plt.figure(figsize=(8, 5))
plt.plot(t, y, label='y(t)', color='blue', linewidth=2)
plt.axvline(t_n, linestyle='--', color='gray', linewidth=1.5)
plt.axhline(y_upper, linestyle='-.', linewidth=1.5, color='gray')
plt.axhline(y_lower, linestyle='-.', linewidth=1.5, color='gray')
plt.text(0.02, 1.05 + 0.02, '±5% предел', va='bottom', fontsize=10)
plt.xlim(0, 2.5)
plt.ylim(0, 1.2)
plt.xlabel('t')
plt.ylabel('y(t)')
plt.title(f'Переходная функция')
plt.grid(True)
plt.tight_layout()
plt.xticks([0, 0.5, 0.9, 1.0, 1.5, 2.0, 2.5])
plt.savefig('images/etalon_response.png', dpi=450)
plt.show()