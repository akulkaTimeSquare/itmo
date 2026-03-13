import numpy as np
from numpy.linalg import norm
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ======================================================
# Вспомогательные функции
# ======================================================

def S(w):
    return np.array([
        [0,     -w[2],  w[1]],
        [w[2],  0,     -w[0]],
        [-w[1], w[0],   0    ]
    ])

# ======================================================
# Параметры винтовой траектории
# ======================================================
R = 15.0      # Радиус цилиндра
PITCH = 5.0   # Шаг винтовой линии

# ======================================================
# Поверхности и градиенты (заданная траектория)
# ======================================================

def phi1(x):
    """Уравнение цилиндра"""
    X, Y, Z = x
    return X**2 + Y**2 - R**2

def phi2(x):
    """Уравнение винтовой линии с коррекцией 'перескока'"""
    X, Y, Z = x
    theta_raw = np.arctan2(Y, X)
    k_revolutions = np.round((Z / PITCH - theta_raw) / (2 * np.pi))
    theta_unwrapped = theta_raw + k_revolutions * 2 * np.pi
    return Z - PITCH * theta_unwrapped

def grad_phi1(x):
    """Градиент уравнения цилиндра"""
    X, Y, Z = x
    return np.array([2*X, 2*Y, 0.0])

def grad_phi2(x):
    """Градиент уравнения винтовой линии"""
    X, Y, Z = x
    den = X**2 + Y**2 + 1e-9
    return np.array([PITCH * Y / den, -PITCH * X / den, 1.0])

# ======================================================
# Параметры системы
# ======================================================

m = 1.2
J = 1.0 * np.eye(3)

Kphi = 1.5
kw = 2.0
V_S_DES = 3.0               # ЖЕЛАЕМАЯ касательная скорость

K = np.diag([4, 4, 4, 3, 3, 3])

nd = np.array([1, 1, 1]) / np.sqrt(3)

# ======================================================
# Численная производная виртуальной скорости
# ======================================================

def vd_func(x):
    g1 = grad_phi1(x)
    g2 = grad_phi2(x)

    # *** КЛЮЧЕВОЕ ИЗМЕНЕНИЕ ***
    # Меняем порядок в векторном произведении, чтобы изменить направление
    # касательного вектора tau и заставить робота двигаться ВВЕРХ.
    tau = np.cross(g2, g1)

    tau_norm = norm(tau)
    tau_hat = tau / (tau_norm + 1e-6) if tau_norm > 0 else tau

    v_att = -Kphi * (phi1(x)*g1 + phi2(x)*g2)
    v_tan = V_S_DES * tau_hat
    return v_att + v_tan

def vd_dot(x, v, h=1e-5):
    return (vd_func(x + h*v) - vd_func(x)) / h

# ======================================================
# Динамика системы (метод пассификации)
# ======================================================

def dynamics(t, state):
    x = state[0:3]
    v = state[3:6]
    n = state[6:9]
    w = state[9:12]

    n = n / norm(n)

    vd = vd_func(x)
    wd = -kw * np.cross(n, nd)

    y = np.hstack([v - vd, w - wd])
    u = -K @ y
    uF = u[0:3]
    uM = u[3:6]

    F = m * (vd_dot(x, v) + uF)
    Mc = uM

    dx = v
    dv = F / m
    dn = S(w) @ n
    dw = np.linalg.inv(J) @ (Mc - np.cross(w, J @ w))

    return np.hstack([dx, dv, dn, dw])

# ======================================================
# Начальные условия
# ======================================================

# Стартуем ниже, чтобы было наглядно видно движение вверх
x0 = np.array([R + 3.0, 2.0, -10.0])
v0 = np.zeros(3)

n0 = np.array([1.0, 1.0, 0.0])
n0 /= norm(n0)

w0 = np.zeros(3)

state0 = np.hstack([x0, v0, n0, w0])

# ======================================================
# Моделирование
# ======================================================

sol = solve_ivp(
    dynamics,
    [0, 40],
    state0,
    max_step=0.01,
    method='RK45'
)

print("Моделирование завершено")

# ======================================================
# График траекторий
# ======================================================

X, Y, Z = sol.y[0], sol.y[1], sol.y[2]

theta = np.linspace(-0.5*np.pi, 3*np.pi, 800)
Xd = R * np.cos(theta)
Yd = R * np.sin(theta)
Zd = PITCH * theta

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

ax.plot(X, Y, Z, label="Фактическая траектория", color='blue', linewidth=2.5)
ax.plot(Xd, Yd, Zd, '--', label="Желаемая траектория (вверх)", color='orange', linewidth=2)
ax.scatter(x0[0], x0[1], x0[2], c='green', s=120, label="Старт", depthshade=True)

cyl_z = np.linspace(min(Zd.min(), Z.min()) - 10, max(Zd.max(), Z.max()) + 10, 20)
cyl_theta = np.linspace(0, 2 * np.pi, 50)
cyl_theta_grid, cyl_z_grid = np.meshgrid(cyl_theta, cyl_z)
cyl_x_grid = R * np.cos(cyl_theta_grid)
cyl_y_grid = R * np.sin(cyl_theta_grid)
ax.plot_surface(cyl_x_grid, cyl_y_grid, cyl_z_grid, alpha=0.1, color='gray')

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title("Слежение за винтовой траекторией")
ax.legend()
plt.tight_layout()
plt.savefig("trajectory_3d.png", dpi=300)
plt.show()



t = sol.t
X, Y, Z = sol.y[0], sol.y[1], sol.y[2]
V = sol.y[3:6].T

# Геометрические ошибки (поверхности)
phi1_err = np.array([phi1(sol.y[0:3, i]) for i in range(len(t))])
phi2_err = np.array([phi2(sol.y[0:3, i]) for i in range(len(t))])

# Полная ошибка положения (евклидова)
pos_error = np.sqrt(phi1_err**2 + phi2_err**2)

# Ошибка скорости
vd_all = np.array([vd_func(sol.y[0:3, i]) for i in range(len(t))])
vel_error = np.linalg.norm(V - vd_all, axis=1)


# ======================================================
# Графики ошибок
# ======================================================

plt.figure()
plt.plot(t, pos_error, linewidth=2)
plt.xlabel("t, c")
plt.ylabel("Ошибка положения")
plt.title("Ошибка слежения e(t) за траекторией")
plt.grid(True)
plt.tight_layout()
plt.savefig("position_error.png", dpi=300)
plt.show()

plt.figure()
plt.plot(t, vel_error, linewidth=2, color='red')
plt.xlabel("t, c")
plt.ylabel("Ошибка скорости")
plt.title("Ошибка слежения по скорости")
plt.grid(True)
plt.tight_layout()
plt.savefig("velocity_error.png", dpi=300)
plt.show()
