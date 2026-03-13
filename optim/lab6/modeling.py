import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.linalg import eig, inv

# --- 1. Исходные данные ---
A = np.array([[3., 1.], [1., 1.]])
B = np.array([[1.], [2.]])
Bf = np.array([[3.], [4.]])
Q = np.array([[3., 0.], [0., 5.]])
n = A.shape[0]

# --- 2. H∞-синтез: Поиск стабилизирующего регулятора ---

# Выбираем gamma > gamma_min. Возьмем с запасом, например, 10.
gamma = 3.5
print(f"Выбранное значение для синтеза: γ = {gamma}\n")

# Решаем уравнение Риккати для выбранного gamma
R_eff = B @ B.T - (1/gamma**2) * Bf @ Bf.T
Hamiltonian = np.block([
    [A, -R_eff],
    [-Q, -A.T]
])

# Находим собственные значения и векторы гамильтоновой матрицы
eigenvalues, eigenvectors = eig(Hamiltonian)

# Отбираем собственные векторы, соответствующие устойчивым полюсам (Re(λ) < 0)
stable_eigenvectors = []
for i, val in enumerate(eigenvalues):
    if np.real(val) < 0:
        stable_eigenvectors.append(eigenvectors[:, i])

# Составляем матрицу U из стабильных собственных векторов
if len(stable_eigenvectors) != n:
    raise ValueError(f"Не удалось найти {n} стабильных собственных векторов. "
                     f"Найдено только {len(stable_eigenvectors)}. Попробуйте другое γ.")

U = np.array(stable_eigenvectors).T
U1 = U[0:n, :]
U2 = U[n:, :]

# Находим матрицу P
P = np.real(U2 @ inv(U1))
print(P)
# Проверяем, что P - положительно определённая
P_eigvals = np.linalg.eigvals(P)
"""print(f"Собственные значения матрицы P: {P_eigvals}")
if np.all(P_eigvals > 0):
    print(">>> Матрица P является положительно определённой.\n")
else:
    print(">>> ВНИМАНИЕ: Матрица P не является положительно определённой.\n")"""


# Вычисляем стабилизирующий регулятор K
K = -B.T @ P
print(f"Рассчитанная матрица регулятора K:\n{K}\n")


# --- 3. Проверка устойчивости и моделирование ---

# Матрица замкнутой системы
A_cl = A + B @ K

# Проверяем устойчивость
cl_eigvals, _ = np.linalg.eig(A_cl)
print(f"Собственные значения замкнутой системы (A + BK): {cl_eigvals}")
if np.all(np.real(cl_eigvals) < 0):
    print(">>> Система УСТОЙЧИВА.\n")
else:
    print(">>> Система НЕУСТОЙЧИВА.\n")

# Начальные условия и возмущение
x0 = np.array([1.0, 0.0])
def disturbance(t):
    return (10 * np.sin(6 * t) + 5 * np.cos(2 * t) +
            4 * np.cos(3 * t) + 3 * np.cos(8 * t))

# Функция динамики системы для решателя ОДУ
def system_dynamics(x, t, A_cl_mat, Bf_vec):
    f_t = disturbance(t)
    dxdt = A_cl_mat @ x + Bf_vec.flatten() * f_t
    return dxdt

# Временной вектор и решение ОДУ
t_sim = np.linspace(0, 15, 1500)
solution = odeint(system_dynamics, x0, t_sim, args=(A_cl, Bf))

x1 = solution[:, 0]
x2 = solution[:, 1]
u = (K @ solution.T).flatten()


# --- 4. Построение графиков ---
plt.figure()
plt.plot(t_sim, x1, label='$x_1(t)$')
plt.plot(t_sim, x2, label='$x_2(t)$')
plt.title('Переменные состояния системы')
plt.ylabel('$x(t)$')
plt.xlabel('t')
plt.grid(True)
plt.legend()
plt.savefig('images/x.png')

plt.figure()
plt.plot(t_sim, u, label='u(t)')
plt.xlabel('t')
plt.ylabel('u(t)')
plt.grid(True)
plt.title('Формируемое управление')
plt.legend()
plt.savefig('images/u.png')

plt.show()
