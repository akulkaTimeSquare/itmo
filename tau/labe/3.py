import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Матрица генератора (блочно-диагональная для частот 1, 2, 3, 7)
G = np.array([[0, 1, 0, 0, 0, 0, 0, 0],
              [-1, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 2, 0, 0, 0, 0],
              [0, 0, -2, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 3, 0, 0],
              [0, 0, 0, 0, -3, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 7],
              [0, 0, 0, 0, 0, 0, -7, 0]])

# Матрицы выходов для разных сигналов
# w = [cos(1t), sin(1t), cos(2t), sin(2t), cos(3t), sin(3t), cos(7t), sin(7t)]
Yg = np.array([[0, 0, 2, 0, 0, 0, 0, 0],
               [0, 0, 0, -3, 0, 0, 0, 0]])  # g = [2*cos(2t), 3*sin(2t)] - исправлена фаза для sin(2t)

Y1 = np.array([[0, 0, 0, 0, 0, 0, 7, 0],
               [0, -3, 0, 0, 0, 0, 0, 0]])  # f1 = [7*cos(7t), 3*sin(1t)] - исправлена фаза для sin(1t)

Y2 = np.array([[0, 0, 0, 0, 0, -5, 0, 0],
               [0, 0, 0, 0, 0, 0, 3, 0]])  # f2 = [5*sin(3t), 3*cos(7t)] - исправлена фаза для sin(3t)

# Единое начальное условие w(0) для генерации всех сигналов
# Состояния: [cos(1t), sin(1t), cos(2t), sin(2t), cos(3t), sin(3t), cos(7t), sin(7t)]
# Начальные условия только 0 или 1:
# - cos(ωt) при t=0: cos(0) = 1, sin(0) = 0 → [1, 0]
# - sin(ωt) при t=0: cos(0) = 0, sin(0) = 1 → [0, 1]
# Для получения нужных амплитуд используем матрицы выходов Y

w0 = np.array([1, 0, 1, 0, 1, 0, 1, 0])  # Единое начальное состояние генератора (только 0 и 1)

def generator_dynamics(x, t):
    return G @ x

# Время моделирования
t = np.linspace(0, 3, 1000)

# Моделирование с единым начальным условием
print("Моделирование внешних воздействий с единым w(0):")

# Решение системы генератора с единым начальным условием
w_t = odeint(generator_dynamics, w0, t)

# Извлечение всех сигналов из единого решения
g_sim = np.array([Yg @ w_t[i] for i in range(len(t))]).T
f1_sim = np.array([Y1 @ w_t[i] for i in range(len(t))]).T
f2_sim = np.array([Y2 @ w_t[i] for i in range(len(t))]).T

plt.figure(figsize=(12, 8))

plt.subplot(2, 3, 1)
plt.plot(t, g_sim[0], 'b-', label='2*cos(2t) симуляция')
plt.plot(t, 2*np.cos(2*t), 'r--', label='2*cos(2t) точное')
plt.title('g: компонента 1')
plt.legend()
plt.grid(True)

plt.subplot(2, 3, 2)
plt.plot(t, g_sim[1], 'b-', label='3*sin(2t) симуляция')
plt.plot(t, 3*np.sin(2*t), 'r--', label='3*sin(2t) точное')
plt.title('g: компонента 2')
plt.legend()
plt.grid(True)

plt.subplot(2, 3, 3)
plt.plot(t, f1_sim[0], 'b-', label='7*cos(7t) симуляция')
plt.plot(t, 7*np.cos(7*t), 'r--', label='7*cos(7t) точное')
plt.title('f1: компонента 1')
plt.legend()
plt.grid(True)

plt.subplot(2, 3, 4)
plt.plot(t, f1_sim[1], 'b-', label='3*sin(1t) симуляция')
plt.plot(t, 3*np.sin(1*t), 'r--', label='3*sin(1t) точное')
plt.title('f1: компонента 2')
plt.legend()
plt.grid(True)

plt.subplot(2, 3, 5)
plt.plot(t, f2_sim[0], 'b-', label='5*sin(3t) симуляция')
plt.plot(t, 5*np.sin(3*t), 'r--', label='5*sin(3t) точное')
plt.title('f2: компонента 1')
plt.legend()
plt.grid(True)

plt.subplot(2, 3, 6)
plt.plot(t, f2_sim[1], 'b-', label='3*cos(7t) симуляция')
plt.plot(t, 3*np.cos(7*t), 'r--', label='3*cos(7t) точное')
plt.title('f2: компонента 2')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('generator_verification.png', dpi=150)
plt.show()

print("\nМатрицы:")
print("G (генератор):")
print(G)
print("\nYg (для g):")
print(Yg)
print("\nY1 (для f1):")
print(Y1)
print("\nY2 (для f2):")
print(Y2)

print("\nЕдиное начальное условие:")
print("w(0) =", w0)
print("Состояния: [cos(1t), sin(1t), cos(2t), sin(2t), cos(3t), sin(3t), cos(7t), sin(7t)]")