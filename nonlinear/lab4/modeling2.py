import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Параметры регулятора
k1 = 5.0
k2 = 5.0

# Моделирвоание
def backstepping_system(t, x):
    x1, x2 = x

    # Виртуальное управление
    a = x1**3 - k1*x1
    z2 = x2 - a

    # Производная a
    da = (3*x1**2 - k1) * (-k1*x1 + z2)

    # Реальное управление
    u = -2*x1 - k2*z2 + da

    # Нелинейная система
    dx1 = -x1**3 + x2
    dx2 = x1 + u

    return [dx1, dx2]

# Начальные условия
x0 = [3, -4]
t_span = (0, 3)
t_eval = np.linspace(*t_span, 2000)

sol = solve_ivp(backstepping_system, t_span, x0, t_eval=t_eval)

x1 = sol.y[0]
x2 = sol.y[1]
t  = sol.t

a = x1**3 - k1*x1
z2 = x2 - a

da = (3*x1**2 -k1) * (-k1*x1 + z2)
u = -2*x1 - k2*z2 + da

# Графики
plt.figure()
plt.plot(t, x1, label=r"$x_1$")
plt.plot(t, x2, label=r"$x_2$")
plt.grid()
plt.legend()
plt.title("Состояния системы")
plt.yticks([-16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 3, 4])
plt.xlabel('t')
plt.ylabel('x(t)')
plt.savefig(
    "images/x_2.png",
    dpi=400,
    bbox_inches="tight",
    pad_inches=0.05
)

plt.figure()
plt.plot(t, z2)
plt.grid()
plt.title(r"Ошибка $z_2$")
plt.xlabel('t')
plt.ylabel('$z_2(t)$')
plt.savefig(
    "images/z2_2.png",
    dpi=400,
    bbox_inches="tight",
    pad_inches=0.05
)

plt.figure()
plt.plot(t, u)
plt.grid()
plt.title("Управление $u$")
plt.xlabel('t')
plt.ylabel('u(t)')
plt.savefig(
    "images/u_2.png",
    dpi=400,
    bbox_inches="tight",
    pad_inches=0.05
)

plt.show()