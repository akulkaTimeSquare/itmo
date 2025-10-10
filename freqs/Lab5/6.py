import numpy as np
import matplotlib.pyplot as plt

def fourier(t, y):
    N = len(t)
    dt = t[1] - t[0]
    m = np.arange(N)
    c = np.sqrt(N) * dt * (-1)**m
    f_numeric = np.fft.fftshift(c * np.fft.fft(y, norm="ortho"))
    return f_numeric

b = 10*np.pi                      # параметр для sinc(bt)
dt = 1e-3
T_true = 200
t = np.arange(-T_true/2, T_true/2, dt)

# ----------------------
# Сигналы
# ----------------------
y2 = np.sinc(b * t / np.pi)

# ----------------------
# Сэмплирование
# ----------------------
Fs = 20                        # частота дискретизации в Гц
dt_sample = 1 / Fs
T = 2
t_sample = np.arange(-T/2, T/2, dt_sample)

y2_sample = np.sinc(b * t_sample / np.pi)

# ----------------------
# Интерполяция через sinc
# ----------------------
def sinc_interp(t, t_sample, y_sample, dt):
    result = np.zeros_like(t)
    for n in range(len(t_sample)):
        result += y_sample[n] * np.sinc((t - t_sample[n])/dt)
    return result

y2_interp = sinc_interp(t, t_sample, y2_sample, dt_sample)


plt.figure(figsize=(14, 8))
font_size = 16
line_width = 6

p = 1
if p == 1:
    plt.plot(t, y2, linewidth=line_width-2, label="Непрерывная функция")
    plt.plot(t, y2_interp, linewidth=line_width+4, linestyle="--", label="Интерполяция")
    plt.stem(t_sample, y2_sample, markerfmt='ro',  linefmt='r-', basefmt=' ', label='Сэмплирование')
    plt.title('Восстановление через теорему Найквиста-Шеннона-Котельникова', fontsize=font_size + 4)
    plt.xlabel('t', fontsize=font_size)
    plt.xlim([-T/2, T/2])
    plt.grid(True)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.legend(fontsize=font_size-1)
    plt.tight_layout()
    plt.savefig(f"images/2_inter_{T}_{Fs}.png", dpi=300)
    plt.show()
elif p == 2:
    plt.plot(t, y2, linewidth=line_width)
    plt.title('Исходная функция', fontsize=font_size + 4)
    plt.xlabel('t', fontsize=font_size)
    plt.xlim([-T/2, T/2])
    plt.grid(True)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.tight_layout()
    plt.savefig(f"images/y2.png", dpi=300)
    plt.show()
elif p == 3:
    V_true = 1/dt
    dnu_true = 1/T_true
    nu_true = np.arange(-V_true/2, V_true/2, dnu_true)
    f2 = fourier(t, y2)

    V_sample = 1/dt_sample
    dnu_sample = 1/T
    nu_sample = np.arange(-V_sample/2, V_sample/2, dnu_sample)
    f_sample = fourier(t_sample, y2_sample)

    f_inter = fourier(t, y2_interp)
    plt.plot(nu_true, f_inter, linewidth=line_width+4, label="Образ интерполяции")
    plt.plot(nu_sample, f_sample, linewidth=line_width, label="Образ сэмплированной")
    plt.plot(nu_true, f2, linewidth=line_width-2, label="Образ непрерывной функции")
    plt.title('Образы Фурье: теорема Найквиста-Шеннона-Котельникова', fontsize=font_size + 4)
    plt.axvline(x=-5, color='red', linestyle='--', linewidth=4)
    plt.axvline(x=5, color='red', linestyle='--', linewidth=4)
    plt.xlabel('v', fontsize=font_size)
    plt.grid(True)
    plt.xlim([-10, 10])
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.legend(fontsize=font_size-1)
    plt.tight_layout()
    plt.savefig(f"images/fourier2_{T}_{Fs}.png", dpi=300)
    plt.show()