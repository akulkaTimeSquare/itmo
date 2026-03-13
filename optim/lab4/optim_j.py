import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import solve_continuous_are, solve_continuous_lyapunov
from scipy.integrate import solve_ivp

# ===============================
# 1. Исходные данные
# ===============================
A = np.array([[0, 1],
              [5, -4]])
b = np.array([[1],
              [1]])
x0 = np.array([1, 0])

Q = np.array([[2, 0],
              [0, 1]])
r = 5
R = np.array([[r]])

# ===============================
# 2. Оптимальное решение LQR
# ===============================
P_opt = solve_continuous_are(A, b, Q, R)
K_opt = np.linalg.inv(R) @ b.T @ P_opt
k1_opt, k2_opt = K_opt[0]
J_opt_val = x0.T @ P_opt @ x0

# ===============================
# 3. Эксперименты
# ===============================
delta = 0.35
experiments = [
    {"name": "k1 + Δ", "k1": k1_opt + delta, "k2": k2_opt, "style": "--", "lw": 1.5},
    {"name": "k1 - Δ", "k1": k1_opt - delta, "k2": k2_opt, "style": "--", "lw": 1.5},
    {"name": "k2 + Δ", "k1": k1_opt, "k2": k2_opt + delta, "style": "--", "lw": 1.5},
    {"name": "k2 - Δ", "k1": k1_opt, "k2": k2_opt - delta, "style": "--", "lw": 1.5},
    {"name": "Оптимальный", "k1": k1_opt, "k2": k2_opt, "style": "k-", "lw": 3},
]

results_data = []

# ===============================
# 4. Подготовка графиков
# ===============================
fig_J, ax_J = plt.subplots()
fig_u, ax_u = plt.subplots()
fig_x, (ax_x1, ax_x2) = plt.subplots(2, 1, sharex=True)

# ===============================
# 5. Основной цикл
# ===============================
for exp in experiments:
    K = np.array([[exp["k1"], exp["k2"]]])
    A_cl = A - b @ K

    # ---- Точное значение J
    Q_aug = Q + K.T @ R @ K
    P_lyap = solve_continuous_lyapunov(A_cl.T, -Q_aug)
    J_val = x0.T @ P_lyap @ x0
    diff_percent = (J_val - J_opt_val) / J_opt_val * 100

    results_data.append({
        "Случай": exp["name"],
        "k1": exp["k1"],
        "k2": exp["k2"],
        "J": J_val,
        "ΔJ (%)": diff_percent
    })

    # ---- Динамика
    def dynamics(t, z):
        x = z[:2].reshape(-1, 1)
        u = -(K @ x)[0, 0]
        dx = A @ x + b * u
        dJ = x.T @ Q @ x + r * u**2
        return [dx[0, 0], dx[1, 0], dJ[0, 0]]

    sol = solve_ivp(
        dynamics,
        [0, 3],
        [x0[0], x0[1], 0],
        t_eval=np.linspace(0, 3, 300)
    )

    t = sol.t
    x1, x2 = sol.y[0], sol.y[1]
    u = -(exp["k1"] * x1 + exp["k2"] * x2)

    # ---- J(t)
    ax_J.plot(t, sol.y[2], exp["style"], lw=exp["lw"],
              label=f"{exp['name']} (J={J_val:.3f})")

    # ---- u(t)
    ax_u.plot(t, u, exp["style"], lw=exp["lw"], label=exp["name"])

    # ---- x(t)
    ax_x1.plot(t, x1, exp["style"], lw=exp["lw"], label=exp["name"])
    ax_x2.plot(t, x2, exp["style"], lw=exp["lw"], label=exp["name"])

# ===============================
# 6. Оформление графиков
# ===============================
ax_J.axhline(J_opt_val, color="black", linestyle="--",
             label=f"J(∞) ≈ {J_opt_val:.2f}")
ax_J.set_title("Накопление функционала J(t)")
ax_J.set_xlabel("t")
ax_J.set_ylabel("J(t)")
ax_J.grid(True)
ax_J.legend()
fig_J.tight_layout()
fig_J.savefig("images/J_comparison.png")

ax_u.set_title("Сравнение управления u(t)")
ax_u.set_xlabel("t")
ax_u.set_ylabel("u(t)")
ax_u.grid(True)
ax_u.legend()
fig_u.tight_layout()
fig_u.savefig("images/u_comparison.png")

ax_x1.set_title("Сравнение состояний")
ax_x1.set_ylabel("x₁(t)")
ax_x2.set_ylabel("x₂(t)")
ax_x2.set_xlabel("t")
ax_x1.grid(True)
ax_x2.grid(True)
ax_x1.legend()
fig_x.tight_layout()
fig_x.savefig("images/x_comparison.png")

plt.show()

# ===============================
# 7. Таблица результатов
# ===============================
df = pd.DataFrame(results_data)
print(df.to_markdown(index=False, float_format="{:.4f}".format))
