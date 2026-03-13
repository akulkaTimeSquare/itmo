import numpy as np
from scipy.linalg import svdvals, eig, inv

# --- 1. Исходные данные и стабилизирующий регулятор (из предыдущих шагов) ---
A = np.array([[3., 1.], [1., 1.]])
B = np.array([[1.], [2.]])
Bf = np.array([[3.], [4.]])
Q = np.array([[3., 0.], [0., 5.]])

# Копируем код для расчета K при gamma = 10, чтобы сделать скрипт самодостаточным
gamma = 3.5
R_eff = B @ B.T - (1/gamma**2) * Bf @ Bf.T
Hamiltonian = np.block([[A, -R_eff], [-Q, -A.T]])
eigenvalues, eigenvectors = eig(Hamiltonian)
stable_eigenvectors = [v for v, l in zip(eigenvectors.T, eigenvalues) if np.real(l) < 0]
U = np.array(stable_eigenvectors).T
P = np.real(U[A.shape[0]:, :] @ inv(U[:A.shape[0], :]))
K = -B.T @ P

# Матрица устойчивой замкнутой системы
A_cl = A + B @ K

# --- 2. Функция для вычисления H∞-нормы ---
def calculate_h_inf_norm(A_matrix, B_matrix, C_matrix, D_matrix, omega_grid):
    """
    Вычисляет H∞-норму для заданной системы (A, B, C, D)
    путем поиска максимального сингулярного значения по сетке частот.
    """
    max_singular_value = 0.0
    identity = np.eye(A_matrix.shape[0])
    
    for w in omega_grid:
        if w == 0: continue # Пропускаем w=0, чтобы избежать деления на 0 при инверсии
        s = 1j * w
        
        # Расчет частотной характеристики G(jω) = C(jωI - A)⁻¹B + D
        try:
            resolvent = inv(s * identity - A_matrix)
            G_jw = C_matrix @ resolvent @ B_matrix + D_matrix
            
            # svdvals возвращает сингулярные значения в убывающем порядке.
            # Берем первое (максимальное).
            current_max_sv = svdvals(G_jw)[0]
            
            if current_max_sv > max_singular_value:
                max_singular_value = current_max_sv
                
        except np.linalg.LinAlgError:
            # Эта ошибка может возникнуть, если частота w совпадает с мнимой частью
            # собственного значения, что для устойчивой системы не должно происходить.
            # Но на всякий случай обработаем.
            print(f"Предупреждение: Ошибка вычисления на частоте w = {w}")
            continue
            
    return max_singular_value

# --- 3. Выполнение расчетов ---

# Создаем достаточно плотную логарифмическую сетку частот для поиска максимума
omega_range = np.logspace(-2, 4, 4000)

# --- Задание 5: H∞-нормы от f до x1 и x2 ---
C1 = np.array([[1., 0.]])
C2 = np.array([[0., 1.]])
D_siso = np.array([[0.]]) # SISO: Single-Input Single-Output

h_inf_norm_fx1 = calculate_h_inf_norm(A_cl, Bf, C1, D_siso, omega_range)
h_inf_norm_fx2 = calculate_h_inf_norm(A_cl, Bf, C2, D_siso, omega_range)

print("--- Результаты вычислений ---")
print("\n[Задание 5]: H∞-нормы для передаточных функций от f до компонент состояния")
print(f"||G_x1_f(s)||∞ = ||C₁ (sI - (A+BK))⁻¹ Bf||∞ ≈ {h_inf_norm_fx1:.4f}")
print(f"||G_x2_f(s)||∞ = ||C₂ (sI - (A+BK))⁻¹ Bf||∞ ≈ {h_inf_norm_fx2:.4f}")


# --- Задание 6: H∞-норма от f до вектора x ---
C_identity = np.eye(A_cl.shape[0])
D_miso = np.zeros((A_cl.shape[0], 1)) # MISO: Multi-Input Single-Output

h_inf_norm_fx = calculate_h_inf_norm(A_cl, Bf, C_identity, D_miso, omega_range)

print("\n[Задание 6]: H∞-норма для передаточной функции от f до вектора состояния")
print(f"||G_x_f(s)||∞  = || (sI - (A+BK))⁻¹ Bf||∞  ≈ {h_inf_norm_fx:.4f}")

