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
dt = 0.01
N = int(T/dt)
t_int = np.linspace(-T/2-dt, T/2+dt, N)

V = 1/dt
dnu = 1/T
nu_int = np.linspace(-V/2-dnu, V/2+dnu, N)

t = np.linspace(-T/2, T/2, 100000)
nu = np.linspace(-V/2, V/2, 100000)


m = np.arange(N)
c = np.sqrt(N) * dt * (-1)**m
f_numeric = np.fft.fftshift(c * np.fft.fft(rect(t_int), norm="ortho"))
rect_numeric = np.fft.ifft(np.fft.ifftshift(f_numeric) / c, norm="ortho")

f_numeric_fft = np.fft.fftshift(np.fft.fft(rect(t_int), norm="ortho"))
rect_numeric_fft = np.fft.ifft(np.fft.ifftshift(f_numeric_fft), norm="ortho")

nu_trapz = np.linspace(-V/2-2*dnu, V/2+2*dnu, N)
t_trapz = np.linspace(-T/2-4*dt, T/2+4*dt, N)
f_numeric_trapz = []
for freq in nu_int:
    integrand = rect(t_int) * np.exp(-2j * np.pi * freq * t_trapz)
    val = np.trapz(integrand, t_trapz)
    f_numeric_trapz.append(val)

rect_numeric_trapz = []
for time in t_int:
    integrand = f_numeric_trapz * np.exp(2j * np.pi * time * nu_trapz)
    val = np.trapz(integrand, nu_trapz)
    rect_numeric_trapz.append(val)

# Построение графиков
plt.figure(figsize=(14, 8))
font_size = 16
line_width = 6

p = 1
if p == 1:
    plt.plot(t, rect(t), color="black", linewidth=line_width+12, label="Аналитическая П(t)")
    plt.plot(t_int, rect_numeric, linewidth=line_width+4, color="r", linestyle="--", label="Восстановление DFT")
    plt.plot(t_int, rect_numeric_fft, linewidth=line_width+6, color="b", linestyle="--", label="Восстановление FFT")
    plt.plot(t_trapz, rect_numeric_trapz, linewidth=line_width+3, color="g", label="Восстановление trapz")
    plt.title('Сравнение методов: восстановление', fontsize=font_size + 4)
    plt.xlabel('t', fontsize=font_size)
    plt.grid(True)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.xlim([-T/2, T/2])
    plt.ylim([-0.15, 1.35])
    plt.legend(fontsize=font_size-1)
    plt.tight_layout()
    plt.savefig(f"images/time_diff_{T}_{dt}.png", dpi=300)
    plt.show()
else:
    plt.plot(nu, fourier_transform(nu), color="black", linewidth=line_width+10, label="Анатилический образ Фурье")
    plt.plot(nu_int, f_numeric, linewidth=line_width+4, color="r", linestyle="-",label="Образ Фурье при DTFT")
    plt.plot(nu_int, f_numeric_fft, linewidth=line_width+4, color="b", linestyle="-",label="Образ Фурье при fft")
    plt.plot(nu_trapz, f_numeric_trapz, linewidth=line_width, color="g", linestyle="--",label="Образ Фурье при trapz")
    plt.title('Сравнение методов: образы', fontsize=font_size + 4)
    plt.xlabel('v', fontsize=font_size)
    plt.grid(True)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.xlim([-V/2, V/2])
#    plt.ylim([-0.35, 1.35])
    plt.legend(fontsize=font_size-1)
    plt.tight_layout()
    #plt.savefig(f"images/freq_diff_w_{T}_{dt}.png", dpi=300)
    plt.show()
