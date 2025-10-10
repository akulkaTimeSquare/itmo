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
dt = 0.15
N = int(T/dt)
#t_int = np.arange(-T/2, T/2+dt, dt)
t_int = np.linspace(-T/2, T/2, N)
#t_int = np.linspace(-T/2-dt, T/2+dt, N)

V = 10
dnu = 0.5
nu_int = np.arange(-V/2, V/2+dnu, dnu)

t = np.linspace(-T/2, T/2, 100000)
nu = np.linspace(-V/2, V/2, 100000)


f_numeric = []
for freq in nu_int:
    integrand = rect(t_int) * np.exp(-2j * np.pi * freq * t_int)
    val = np.trapz(integrand, t_int)
    f_numeric.append(val)

rect_numeric = []
for time in t_int:
    integrand = f_numeric * np.exp(2j * np.pi * time * nu_int)
    val = np.trapz(integrand, nu_int)
    rect_numeric.append(val)


# Построение графиков
plt.figure(figsize=(14, 8))
font_size = 16
line_width = 6

p = 1
if p == 1:
    plt.plot(t, rect(t), linewidth=line_width, label="Аналитическая П(t)")
    plt.plot(t_int, rect_numeric, linewidth=line_width, color="r", label="Восстановленная функция")
    plt.title('Восстановление при численном интегрировании', fontsize=font_size + 4)
    plt.xlabel('t', fontsize=font_size)
    plt.grid(True)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.xlim([-T/2, T/2])
    plt.ylim([-0.35, 1.45])
    plt.legend(fontsize=font_size-1)
    plt.tight_layout()
    plt.savefig(f"images/time_{T}_{dt}_{V}_{dnu}.png", dpi=300)
    plt.show()
else:
    plt.plot(nu, fourier_transform(nu), linewidth=line_width, label="Анатилический образ Фурье")
    plt.plot(nu_int, f_numeric, linewidth=line_width, color="r", linestyle=":", label="Образ при численном интегрировании")
    plt.title('Образы Фурье при численном интегрировании', fontsize=font_size + 4)
    plt.xlabel('v', fontsize=font_size)
    plt.grid(True)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.xlim([-V/2, V/2])
#    plt.ylim([-0.35, 1.35])
    plt.legend(fontsize=font_size-1)
    plt.tight_layout()
    plt.savefig(f"images/freq_{T}_{dt}_{V}_{dnu}.png", dpi=300)
    plt.show()


