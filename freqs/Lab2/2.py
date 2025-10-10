import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy.integrate import trapezoid

# 1. Загрузка аудиофайла
Fs, y = wavfile.read('ac.mp3')  # замените на имя вашего файла

# 2. Преобразование в моно, если нужно
if y.ndim > 1:
    y = y[:, 0]  # используем только один канал

# 3. Создание временного вектора
N = len(y)
t = np.arange(N) / Fs

# 4. Построение графика f(t)
plt.figure()
plt.plot(t, y)
plt.title('Сигнал во времени')
plt.xlabel('Время [с]')
plt.ylabel('Амплитуда')
plt.grid(True)
plt.show()

# 5. Численное преобразование Фурье (без fft)
V = 1000               # максимальная частота
dv = 1                 # шаг по частоте
v = np.arange(-V, V+1, dv)  # диапазон частот
Y = np.zeros_like(v, dtype=complex)

# 6. Численное интегрирование методом трапеций
for k in range(len(v)):
    integrand = y * np.exp(-1j * 2 * np.pi * v[k] * t)
    Y[k] = trapezoid(integrand, t)

# 7. Построение |F(ν)|
plt.figure()
plt.plot(v, np.abs(Y))
plt.title('Амплитудный спектр |F(ν)|')
plt.xlabel('Частота [Гц]')
plt.ylabel('Амплитуда')
plt.grid(True)
plt.show()