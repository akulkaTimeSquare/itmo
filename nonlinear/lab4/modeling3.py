import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import sympy as sp

# Символьное определение
x1, x2, x3, x4, t = sp.symbols('x1 x2 x3 x4 t')
k1, k2, k3, k4 = sp.symbols('k1 k2 k3 k4', positive=True, real=True)

# Виртуальные управления
a1 = sp.cos(x1) + k1 * x1
z2 = x2 - a1

x1_dot = -k1 * x1 - z2
x1_dot_orig = sp.cos(x1) - x2
x2_dot_orig = x1 + x3
x3_dot_orig = x1 * x3 + (2 - sp.sin(x3)) * x4
x4_dot_orig = x2 * x3 + 2 * sp.Symbol('u')

# Производные
da1 = sp.diff(a1, x1) * x1_dot_orig
da1_simpl = sp.simplify(da1)
print("da1 =", da1_simpl)

a2 = da1_simpl - k2 * z2
da2 = sp.diff(a2, x1) * x1_dot_orig + sp.diff(a2, x2) * x2_dot_orig
da2_simpl = sp.simplify(da2)
print("da2 =", da2_simpl)

z3 = x3 - a2
a3_expr = (-x1 * x3 + da2_simpl - z2 - k3 * z3) / (2 - sp.sin(x3))
a3_simpl = sp.simplify(a3_expr)
print("a3 =", a3_simpl)

da3 = (sp.diff(a3_simpl, x1) * x1_dot_orig +
       sp.diff(a3_simpl, x2) * x2_dot_orig +
       sp.diff(a3_simpl, x3) * x3_dot_orig)
da3_simpl = sp.simplify(da3)
print("da3 (упрощённая) = ... (длинная, но точная)")

# Управление
z4 = x4 - a3_simpl
u_expr = (-x2 * x3 + da3_simpl - z3 - k4 * z4) / 2
u_simpl = sp.simplify(u_expr)
print("u = ... (точное выражение)")

# Параметры и функции
subs_dict = {k1: 2.0, k2: 3.0, k3: 4.0, k4: 5.0}

a1_func   = sp.lambdify((x1,), a1.subs(subs_dict), 'numpy')
da1_func  = sp.lambdify((x1, x2), da1_simpl.subs(subs_dict), 'numpy')
a2_func   = sp.lambdify((x1, x2), a2.subs(subs_dict), 'numpy')
da2_func  = sp.lambdify((x1, x2, x3), da2_simpl.subs(subs_dict), 'numpy')
a3_func   = sp.lambdify((x1, x2, x3), a3_simpl.subs(subs_dict), 'numpy')
da3_func  = sp.lambdify((x1, x2, x3, x4), da3_simpl.subs(subs_dict), 'numpy')
u_func    = sp.lambdify((x1, x2, x3, x4), u_simpl.subs(subs_dict), 'numpy')

# Динамика
def system_dynamics(t, X):
    x1v, x2v, x3v, x4v = X

    try:
        u_val = u_func(x1v, x2v, x3v, x4v)
    except Exception as e:
        print(f"Ошибка в u при t={t:.3f}, x={X}: {e}")
        u_val = 0.0

    dx1 = np.cos(x1v) - x2v
    dx2 = x1v + x3v
    dx3 = x1v * x3v + (2 - np.sin(x3v)) * x4v
    dx4 = x2v * x3v + 2 * u_val

    return [dx1, dx2, dx3, dx4]

x0 = [2, -1, 1, -2]

t_span = (0, 4)
t_eval = np.linspace(*t_span, 1000)

sol = solve_ivp(
    system_dynamics, t_span, x0,
    t_eval=t_eval,
    method='RK45',
    rtol=1e-8,
    atol=1e-10,
    max_step=0.01
)

if not sol.success:
    print("Интегрирование не удалось:", sol.message)
else:
    x1, x2, x3, x4 = sol.y

    a1_vals = a1_func(x1)
    z2_vals = x2 - a1_vals

    da1_vals = da1_func(x1, x2)
    a2_vals = a2_func(x1, x2)
    z3_vals = x3 - a2_vals

    a3_vals = a3_func(x1, x2, x3)
    z4_vals = x4 - a3_vals

    u_vals = u_func(x1, x2, x3, x4)
    u_vals = np.clip(u_vals, -50, 50)


    # Графики
    plt.figure()
    plt.plot(sol.t, x1, label=r'$x_1$')
    plt.plot(sol.t, x2, label=r'$x_2$')
    plt.plot(sol.t, x3, label=r'$x_3$')
    plt.plot(sol.t, x4, label=r'$x_4$')
    plt.grid(True)
    plt.legend()
    plt.title('Состояния системы')
    plt.xlabel('t')
    plt.ylabel('x(t)')
    plt.yticks([-8, -6, -4, -2, -1, 0, 1, 2, 4, 6, 8])
    plt.savefig(
        "images/x_3.png",
        dpi=400,
        bbox_inches="tight",
        pad_inches=0.05
    )

    plt.figure()
    plt.plot(sol.t, z2_vals, label=r'$z_2$')
    plt.plot(sol.t, z3_vals, label=r'$z_3$')
    plt.plot(sol.t, z4_vals, label=r'$z_4$')
    plt.grid(True)
    plt.legend()
    plt.title('Ошибки управления $z_i$')
    plt.xlabel('t')
    plt.ylabel('$z_i(t)$')
    plt.savefig(
        "images/z_3.png",
        dpi=400,
        bbox_inches="tight",
        pad_inches=0.05
    )


    plt.figure()
    plt.plot(sol.t, u_vals)
    plt.grid(True)
    plt.legend()
    plt.title('Управление u')
    plt.xlabel('t')
    plt.ylabel('$u(t)$')
    plt.savefig(
        "images/u_3.png",
        dpi=400,
        bbox_inches="tight",
        pad_inches=0.05
    )

    plt.show()