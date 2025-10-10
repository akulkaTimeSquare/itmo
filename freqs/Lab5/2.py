import numpy as np
import matplotlib.pyplot as plt

# Прямоугольная функция П(t)
def rect(t):
    return np.where(np.abs(t) <= 0.5, 1, 0)

# Фурье-образ (синк-функция)
def fourier_transform(nu):
    return np.sinc(nu)  # np.sinc(x) = sin(pi x)/(pi x)

# Диапазоны значений

T = 5
dt = 0.05
N = int(T/dt)
t_int = np.linspace(-T/2, T/2, N)

V = 1/dt
dnu = 1/T
nu_int = np.linspace(-V/2-dnu, V/2, N)

t = np.linspace(-T/2, T/2, 100000)
nu = np.linspace(-V/2, V/2, 100000)


f_numeric = np.fft.fftshift(np.fft.fft(rect(t_int), norm="ortho"))

rect_numeric = np.fft.ifft(np.fft.ifftshift(f_numeric), norm="ortho")


# Построение графиков
plt.figure(figsize=(14, 8))
font_size = 16
line_width = 6

p = 1
if p == 1:
    plt.plot(t, rect(t), linewidth=line_width, label="Аналитическая П(t)")
    plt.plot(t_int, rect_numeric, linewidth=line_width, color="r", linestyle="-", label="Восстановленная функция")
    plt.title('Восстановление при БПФ', fontsize=font_size + 4)
    plt.xlabel('t', fontsize=font_size)
    plt.grid(True)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.xlim([-T/2, T/2])
#    plt.ylim([-1.3, 1.3])
    plt.legend(fontsize=font_size-1)
    plt.tight_layout()
    plt.savefig(f"images/time_fft_{T}_{dt}.png", dpi=300)
    plt.show()
else:
    plt.plot(nu_int, f_numeric, linewidth=line_width, color="r", linestyle="-",label="Образ Фурье при БПФ")
    plt.plot(nu, fourier_transform(nu), linewidth=line_width+1, label="Анатилический образ Фурье")
    plt.title('Образы Фурье при БПФ', fontsize=font_size + 4)
    plt.xlabel('v', fontsize=font_size)
    plt.grid(True)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.xlim([-V/2, V/2])
#    plt.ylim([-0.35, 1.35])
    plt.legend(fontsize=font_size-1)
    plt.tight_layout()
    plt.savefig(f"images/freq_fft_{T}_{dt}.png", dpi=300)
    plt.show()


