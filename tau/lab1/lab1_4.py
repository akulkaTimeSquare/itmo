# Выполню расчёт и построю графики управления u(t) и эволюции состояния x(t).
import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Данные из задачи
A = np.array([[7, -6, 9],
              [6, -5, 6],
              [-6, 3, -8]], dtype=np.float64)
B = np.array([[2], [1], [-1]], dtype=np.float64)
x1 = np.array([[-2], [1], [-1]], dtype=np.float64)
t1 = 3.0

def integrand(tau):
    expA = expm(A * tau)
    expAT = expm(A.T * tau)
    return expA @ (B @ B.T) @ expAT

# Интегрирование поэлементно
Wc = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        f = lambda tau: integrand(tau)[i, j]
        Wc[i, j], _ = quad(f, 0, 3)

pinvWc = np.linalg.pinv(Wc)

print("Wc:\n", Wc)
print("Pseudoinverse Wc:\n", pinvWc)

# функция управления u(t)
def u_of_t(t):
    # scalar: B^T * exp(A^T*(t1-t)) * P_inv * x1
    E = expm(A.T * (t1 - t))
    return float((B.T @ E @ pinvWc @ x1).squeeze())

# рассчитаем u на множестве точек
ts = np.linspace(0, t1, 500)
us = np.array([u_of_t(t) for t in ts])

# интегрируем систему x' = A x + B u(t), с начальным x(0)=0
def rhs(t, x):
    # x приходят как вектор (3,)
    ut = u_of_t(t)
    return (A @ x.reshape(-1,1) + B * ut).flatten()

x0 = np.zeros(3)
sol = solve_ivp(rhs, (0, t1), x0, t_eval=ts, rtol=1e-8, atol=1e-10)

# построим график u(t)
plt.figure(figsize=(8,4))
plt.plot(ts, us)
plt.xlabel('t')
plt.ylabel('u(t)')
plt.title('Управление $u(t)=B^T e^{A^T (t_1 - t)} (P(t_1))^{+} x_1$')
plt.grid(True)
plt.tight_layout()
plt.savefig('images/lab1_4_1.png')
plt.show()

# график компонент состояния x(t)
plt.figure(figsize=(8, 6))
plt.plot(sol.t, sol.y[0, :], label='x1(t)')
plt.plot(sol.t, sol.y[1, :], label='x2(t)')
plt.plot(sol.t, sol.y[2, :], label='x3(t)')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('Состояние системы x(t) при x(0)=0 и управлении $u(t)=B^T e^{A^T (t_1 - t)} (P(t_1))^{+} x_1$')
plt.legend()
current_yticks = plt.yticks()[0]
additional_ticks = [-2, -1, 1]
updated_yticks = sorted(set(current_yticks) | set(additional_ticks))
plt.yticks(updated_yticks)
plt.grid(True)
plt.tight_layout()
plt.savefig('images/lab1_4_2.png')
plt.show()