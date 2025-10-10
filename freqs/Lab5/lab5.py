import numpy as np
import matplotlib.pyplot as plt


t = np.linspace(-1, 1, 10**2)
y = np.array([])
for i in t:
    y = np.append(y, [1]) if abs(i) <= 1 / 2 else np.append(y, [0])

v = np.linspace(-10, 10, 10**2)
fourier = np.array([])
for i in v:
    fourier = np.append(fourier, [np.sin(np.pi*i)/(np.pi*i)])

num_fourier = np.array([])
for i in v:
    num_fourier = np.append(num_fourier, [np.trapz(y*np.exp(-2*np.pi*1j*i*t), t)])

num_y = np.array([])
for i in t:
    num_y = np.append(num_y, [np.trapz(num_fourier*np.exp(2*np.pi*1j*i*v), v)])

num_fourier_f = np.fft.fftshift(np.fft.fft(y, norm="ortho"))
num_y_f = np.fft.ifft(np.fft.ifftshift(num_fourier_f), norm="ortho")

"""
figure, [ax1, ax2] = plt.subplots(2, 1)
ax1.plot(t, y)
ax1.grid(True)
ax2.plot(v, fourier)
ax2.grid(True)
figure.tight_layout()
figure.savefig('images/1.png')
"""

"""
plt.plot(v, fourier, color='b', label="Истинный Фурье-образ")
plt.plot(v, num_fourier, color='r', linestyle='dashed', linewidth=3, label="Численно вычисленный")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('images/2.png')
"""

"""
plt.plot(t, y, color='b', label="Истинная функция")
plt.plot(t, num_y, color='r', linestyle='dashed', linewidth=3, label="Численно вычисленная")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('images/3.png')
"""

"""
plt.plot(v, fourier, color='b', label="Истинный Фурье-образ")
plt.plot(v, num_fourier_f, color='r', linestyle='dashed', linewidth=3, label="Численно вычисленная (fft)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('images/4.png')
plt.show()
"""


plt.plot(t, y, color='b', label="Истинная функция")
plt.plot(t, num_y_f, color='r', linestyle='dashed', linewidth=3, label="Численно вычисленная (ifft)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('images/5.png')
plt.show()
