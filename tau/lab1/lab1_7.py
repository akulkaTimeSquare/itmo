import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy.integrate import ode

A = np.array([
    [-21, -38, 6],
    [8, 13, -4],
    [-6, -14, -1]
])
T = np.array([
    [2, 3, 2],
    [-1, -1, -1],
    [1, 1, 0]
])

C = np.array([[7, 14, 0]])

print("C@T:", C@T)

# Формируем матрицу наблюдаемости
V = np.vstack([
    C,
    C @ A,
    C @ np.linalg.matrix_power(A, 2)
])


print("Матрица наблюдаемости V:")
print(V)

# Вычисляем ранг матрицы наблюдаемости
rank = np.linalg.matrix_rank(V)
print("Ранг матрицы наблюдаемости V:", rank)

# Вычисляем собственные числа матрицы A
eigenvalues = np.linalg.eigvals(A)
print("Собственные числа матрицы A:", eigenvalues)

H1 = np.vstack([A - eigenvalues[0] * np.eye(3), C])
print("H1:", H1)
H2 = np.vstack([A - eigenvalues[1] * np.eye(3), C])
print("H2:", H2)
H3 = np.vstack([A - eigenvalues[2] * np.eye(3), C])
print("H3:", H3)

print("np.linalg.matrix_rank(H1):", np.linalg.matrix_rank(H1))
print("np.linalg.matrix_rank(H2):", np.linalg.matrix_rank(H2))
print("np.linalg.matrix_rank(H3):", np.linalg.matrix_rank(H3))

def integrand(tau):
    expA = expm(A * tau)
    expAT = expm(A.T * tau)
    return expAT @ (C.T @ C) @ expA

Q = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        f = lambda tau: integrand(tau)[i, j]
        Q[i, j], _ = quad(f, 0, 3)

print("Q:", Q)

eigenvalues = np.linalg.eigvals(Q)
print("eigenvalues:", eigenvalues)


pinvQ = np.linalg.pinv(Q)
print("Q@invQ:", Q@pinvQ@Q)

def y(t):
    return np.array([3*np.exp(-5*t) * np.cos(2*t) - 1*np.exp(-5*t) * np.sin(2*t)])

def integrandx0(tau):
    expAT = expm(A.T * tau)
    return expAT @ C.T @ y(tau)

intx0 = np.zeros((3, 1))
for i in range(3):
    f = lambda tau: integrandx0(tau)[i]
    intx0[i], _ = quad(f, 0, 3)

print("intx0:", intx0)

x0 = pinvQ @ intx0
print("x0:", x0)
v = np.array([[2], [-1], [1]])

# Функция для моделирования системы
def system_dynamics(t, x):
    """Динамика системы dx/dt = Ax"""
    return A @ x

# Построение графиков
def plot_graphs():
    # Временной интервал
    t_start, t_end = 0, 3
    dt = 0.001
    t = np.arange(t_start, t_end, dt)
    
    # Значения alpha для исследования
    alphas = [0, 1, 5]
    
    # Сохраняем результаты для всех alpha
    all_results = {}
    
    for alpha in alphas:
        # Вычисляем начальное условие для данного alpha
        x0_alpha = pinvQ @ intx0 + alpha * v
        
        # Решаем систему dx/dt = Ax с начальным условием x0_alpha
        solver = ode(system_dynamics)
        solver.set_integrator('dopri5')
        solver.set_initial_value(x0_alpha.flatten(), t_start)
        
        x_system = []
        t_system = []
        
        while solver.successful() and solver.t < t_end:
            solver.integrate(solver.t + dt)
            x_system.append(solver.y)
            t_system.append(solver.t)
        
        x_system = np.array(x_system)
        t_system = np.array(t_system)
        
        # Сохраняем результаты
        all_results[alpha] = {
            'x_system': x_system,
            't_system': t_system,
            'x0': x0_alpha
        }
        
        # График 1: Состояние системы x(t) для каждого alpha
        plt.figure(figsize=(10, 6))
        plt.plot(t_system, x_system[:, 0], 'b-', linewidth=2, label='x₁(t)')
        plt.plot(t_system, x_system[:, 1], 'r-', linewidth=2, label='x₂(t)')
        plt.plot(t_system, x_system[:, 2], 'g-', linewidth=2, label='x₃(t)')
        plt.xlabel('t')
        plt.ylabel('x(t)')
        if alpha == 0:
            plt.title('Состояние системы x(t) при начальном состоянии x(0)')
        else:
            plt.title(f'Состояние системы x(t) при начальном состоянии x(0) + {alpha}v')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'images/lab1_7_state{alpha}.png', dpi=300, bbox_inches='tight')
        plt.show()

        # График 2: Сравнение истинного y(t) с выходом системы Cx(t) для каждого alpha
        plt.figure(figsize=(10, 6))
        
        # Истинное значение y(t)
        y_true = [y(t_val)[0] for t_val in t_system]
        
        # Выход системы Cx(t)
        y_system = [C @ x for x in x_system]
        y_system = [y_val[0] for y_val in y_system]

        plt.plot(t_system, y_true, 'b-', linewidth=2, label='f(t) - истинное')
        plt.plot(t_system, y_system, 'r--', linewidth=2, label='y = Cx(t) - выход системы')
        plt.xlabel('t')
        plt.ylabel('y(t)')
        if alpha == 0:
            plt.title('Сравнение y(t) и f(t) при начальном состоянии x(0)')
        else:
            plt.title(f'Сравнение y(t) и f(t) при начальном состоянии x(0) + {alpha}v')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'images/lab1_7_comparison{alpha}.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    # График 3: Сравнение компонент x₁(t) для всех alpha
    plt.figure(figsize=(8, 4))
    
    colors = ['blue', 'red', 'green']
    styles = ['-', '-', '-']
    
    for i, alpha in enumerate(alphas):
        t_system = all_results[alpha]['t_system']
        x_system = all_results[alpha]['x_system']
        
        if alpha == 0:
            label = 'x₁(t) при x(0)'
        else:
            label = f'x₁(t) при x(0) + {alpha}v'
        
        plt.plot(t_system, x_system[:, 0], color=colors[i], linestyle=styles[i], 
                linewidth=2, label=label)
    
    plt.xlabel('t')
    plt.ylabel('x₁(t)')
    plt.title('Сравнение компоненты x₁(t) при различных начальных условиях')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('images/lab1_7_x1_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # График 4: Сравнение компонент x₂(t) для всех alpha
    plt.figure(figsize=(8, 4))
    
    for i, alpha in enumerate(alphas):
        t_system = all_results[alpha]['t_system']
        x_system = all_results[alpha]['x_system']
        
        if alpha == 0:
            label = 'x₂(t) при x(0)'
        else:
            label = f'x₂(t) при x(0) + {alpha}v'
        
        plt.plot(t_system, x_system[:, 1], color=colors[i], linestyle=styles[i], 
                linewidth=2, label=label)
    
    plt.xlabel('t')
    plt.ylabel('x₂(t)')
    plt.title('Сравнение компоненты x₂(t) при различных начальных условиях')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('images/lab1_7_x2_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # График 5: Сравнение компонент x₃(t) для всех alpha
    plt.figure(figsize=(8, 4))
    
    for i, alpha in enumerate(alphas):
        t_system = all_results[alpha]['t_system']
        x_system = all_results[alpha]['x_system']
        
        if alpha == 0:
            label = 'x₃(t) при x(0)'
        else:
            label = f'x₃(t) при x(0) + {alpha}v'
        
        plt.plot(t_system, x_system[:, 2], color=colors[i], linestyle=styles[i], 
                linewidth=2, label=label)
    
    plt.xlabel('t')
    plt.ylabel('x₃(t)')
    plt.title('Сравнение компоненты x₃(t) при различных начальных условиях')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('images/lab1_7_x3_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # График 6: Сравнение выходов всех систем на одном графике
    plt.figure(figsize=(8, 4))
    
    for i, alpha in enumerate(alphas):
        t_system = all_results[alpha]['t_system']
        x_system = all_results[alpha]['x_system']
        
        # Выход системы Cx(t)
        y_system = [C @ x for x in x_system]
        y_system = [y_val[0] for y_val in y_system]
        
        if alpha == 0:
            label = 'y = Cx(t) при x(0)'
        else:
            label = f'y = Cx(t) при x(0) + {alpha}v'
        
        plt.plot(t_system, y_system, color=colors[i], linestyle=styles[i], 
                linewidth=2, label=label)
    
    # Истинное значение y(t) (одинаковое для всех)
    y_true = [y(t_val)[0] for t_val in all_results[0]['t_system']]
    plt.plot(all_results[0]['t_system'], y_true, 'k-', linewidth=3, 
            label='f(t) - истинное', alpha=0.7)
    
    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.title('Сравнение выходов систем при различных начальных условиях')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('images/lab1_7_all_outputs.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Выводим начальные условия для всех alpha
    print("\nНачальные условия для различных alpha:")
    for alpha in alphas:
        x0_alpha = all_results[alpha]['x0']
        print(f"alpha = {alpha}: x(0) = [{x0_alpha[0,0]:.3f}, {x0_alpha[1,0]:.3f}, {x0_alpha[2,0]:.3f}]ᵀ")

# Запуск построения графиков
if __name__ == "__main__":
    plot_graphs()

