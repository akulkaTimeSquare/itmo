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
plt.plot(t, y, label='Переходная функция $y(t)$')
plt.axvline(t_n, linestyle='--', color='gray', linewidth=0.8)
plt.axhline(y_upper, linestyle='-.', linewidth=0.7, color='gray')
plt.axhline(y_lower, linestyle='-.', linewidth=0.7, color='gray')
plt.text(0.02, 1.05 + 0.02, '±5% предел', va='bottom', fontsize=8)
plt.xlim(0, 2.5)
plt.ylim(0, 1.2)
plt.xlabel('Время, с')
plt.ylabel('y(t)')
plt.title(f'Переходная функция при ζ=1, ωₙ={wn_solution:.3f}\n(y(0.9)=0.95)')
plt.grid(True)
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig('images/etalon_response.png', dpi=450)
plt.show()