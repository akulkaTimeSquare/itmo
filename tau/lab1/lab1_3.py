# Выполню расчёт и построю графики управления u(t) и эволюции состояния x(t).
import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd

# Данные из задачи
A = np.array([[7, -6, 9],
              [6, -5, 6],
              [-6, 3, -8]], dtype=float)
B = np.array([[-2], [-1], [2]], dtype=float)
x1 = np.array([[-5], [-3], [3]], dtype=float)
P_given = np.array([[0.3688436, 0.20346076, -0.37345989],
                    [0.20346076, 0.13461396, -0.17461489],
                    [-0.37345989, -0.17461489, 0.46461429]], dtype=float)
t1 = 3.0

# обратная матрица P
P_inv = np.linalg.inv(P_given)

# функция управления u(t)
def u_of_t(t):
    # scalar: B^T * exp(A^T*(t1-t)) * P_inv * x1
    E = expm(A.T * (t1 - t))
    return float((B.T @ E @ P_inv @ x1).squeeze())

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
plt.title('Управление $u(t)=B^T e^{A^T (t_1 - t)} (P(t_1))^{-1} x_1$')
plt.grid(True)
plt.tight_layout()
plt.savefig('images/lab1_3_1.png')
plt.show()

# график компонент состояния x(t)
plt.figure(figsize=(8,4))
plt.plot(sol.t, sol.y[0,:], label='x1(t)')
plt.plot(sol.t, sol.y[1,:], label='x2(t)')
plt.plot(sol.t, sol.y[2,:], label='x3(t)')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('Состояние системы x(t) при x(0)=0 и управлении $u(t)=B^T e^{A^T (t_1 - t)} (P(t_1))^{-1} x_1$')
plt.yticks(np.arange(-5, 5, 1))
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('images/lab1_3_2.png')
plt.show()

# проверим значение x(t1)
x_t1 = sol.y[:, -1].reshape(-1,1)
check = np.hstack([x_t1, x1])
df_check = pd.DataFrame(check, index=['x1','x2','x3'], columns=['x(t1) numeric','x1 target'])
df_P = pd.DataFrame(P_given, index=['r1','r2','r3'], columns=['c1','c2','c3'])
df_Pinv = pd.DataFrame(P_inv, index=['r1','r2','r3'], columns=['c1','c2','c3'])
