import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

def get_eigenvalues(A):
    return np.linalg.eig(A)[0]

def compute_transfer_matrix(A, B, C):
    s = sp.symbols('s')
    A_sp = sp.Matrix(A.tolist())
    B_sp = sp.Matrix(B.tolist())
    C_sp = sp.Matrix(C.tolist())
    I = sp.eye(A_sp.shape[0])
    W = C_sp * (s*I - A_sp).inv() * B_sp
    return sp.simplify(W)

def compute_det_and_zero_pole(W):
    s = sp.symbols('s')
    detW = sp.simplify(sp.factor(W.det()))
    num, den = sp.fraction(detW)
    num = sp.expand(num)
    den = sp.expand(den)
    num_poly = sp.Poly(num, s)
    den_poly = sp.Poly(den, s)
    zeros_exact = sp.roots(num_poly)
    poles_exact = sp.roots(den_poly)
    zeros_numeric = [complex(z) for z in sp.nroots(num_poly)] if num_poly.degree() > 0 else []
    poles_numeric = [complex(p) for p in sp.nroots(den_poly)] if den_poly.degree() > 0 else []
    return detW, zeros_exact, poles_exact, zeros_numeric, poles_numeric

def controllability_matrix(A, B):
    A_sp = sp.Matrix(A.tolist())
    B_sp = sp.Matrix(B.tolist())
    n = A_sp.shape[0]
    blocks = [B_sp]
    Ak = sp.eye(n)
    for _ in range(1, n):
        Ak = Ak * A_sp
        blocks.append(Ak * B_sp)
    return sp.Matrix.hstack(*blocks)

def observability_matrix(A, C):
    A_sp = sp.Matrix(A.tolist())
    C_sp = sp.Matrix(C.tolist())
    n = A_sp.shape[0]
    rows = [C_sp]
    Ak = sp.eye(n)
    for _ in range(1, n):
        Ak = Ak * A_sp
        rows.append(C_sp * Ak)
    return sp.Matrix.vstack(*rows)

def output_controllability_matrix(A, B, C):
    A_sp = sp.Matrix(A.tolist())
    B_sp = sp.Matrix(B.tolist())
    C_sp = sp.Matrix(C.tolist())
    n = A_sp.shape[0]
    blocks = [C_sp * B_sp]
    Ak = sp.eye(n)
    for _ in range(1, n):
        Ak = Ak * A_sp
        blocks.append(C_sp * Ak * B_sp)
    return sp.Matrix.hstack(*blocks)

def pbh_controllable(A, B):
    A_sp = sp.Matrix(A.tolist())
    B_sp = sp.Matrix(B.tolist())
    n = A_sp.shape[0]
    I = sp.eye(n)
    eigs = list(A_sp.eigenvals().keys())
    uncontrollable = []
    for lam in eigs:
        M = (lam*I - A_sp).row_join(B_sp)
        if M.rank() < n:
            uncontrollable.append(lam)
    return len(uncontrollable) == 0, uncontrollable

def pbh_observable(A, C):
    return pbh_controllable(sp.Matrix(A.tolist()).T, sp.Matrix(C.tolist()).T)

def is_stabilizable(A, B):
    ok, uncontrollable = pbh_controllable(A, B)
    if ok:
        return True, []
    unstable_uncontrollable = []
    for lam in uncontrollable:
        if sp.re(complex(lam.evalf())) >= 0:
            unstable_uncontrollable.append(lam)
    return len(unstable_uncontrollable) == 0, unstable_uncontrollable

def is_detectable(A, C):
    ok, unobservable = pbh_observable(A, C)
    if ok:
        return True, []
    unstable_unobservable = []
    for lam in unobservable:
        if sp.re(complex(lam.evalf())) >= 0:
            unstable_unobservable.append(lam)
    return len(unstable_unobservable) == 0, unstable_unobservable

def symbolic_impulse_response(A, B, C):
    t = sp.symbols('t', real=True, positive=True)
    A_sp = sp.Matrix(A.tolist())
    B_sp = sp.Matrix(B.tolist())
    C_sp = sp.Matrix(C.tolist())
    g = sp.simplify(C_sp * sp.exp(A_sp * t) * B_sp)
    return t, g

def symbolic_step_response(A, B, C):
    t = sp.symbols('t', real=True, positive=True)
    A_sp = sp.Matrix(A.tolist())
    B_sp = sp.Matrix(B.tolist())
    C_sp = sp.Matrix(C.tolist())
    I = sp.eye(A_sp.shape[0])
    # h(t) = C A^{-1} (e^{At} - I) B, valid if A is invertible
    h = sp.simplify(C_sp * A_sp.inv() * (sp.exp(A_sp * t) - I) * B_sp)
    return t, h

def evaluate_matrix_function_t(mat_func_sym, t_sym, t_values):
    f = sp.lambdify(t_sym, mat_func_sym, modules=['numpy'])
    values = [np.array(f(tv), dtype=complex) for tv in t_values]
    return np.stack(values, axis=-1)  # shape: (p, m, T)

def frequency_response(A, B, C, w_values):
    A_np = np.array(A, dtype=complex)
    B_np = np.array(B, dtype=complex)
    C_np = np.array(C, dtype=complex)
    n = A_np.shape[0]
    p = C_np.shape[0]
    m = B_np.shape[1]
    Wjw = np.zeros((p, m, len(w_values)), dtype=complex)
    I = np.eye(n, dtype=complex)
    for k, w in enumerate(w_values):
        M = 1j*w*I - A_np
        Winv = np.linalg.inv(M)
        Wjw[:, :, k] = C_np @ Winv @ B_np
    return Wjw

A = np.array([[0, 1], [-1, 1]])
print(get_eigenvalues(A))

B = np.array([[1, 2], [1, 0]])
C = np.array([[1, -2], [0, 3]])

W = compute_transfer_matrix(A, B, C)
print("W(s) =")
sp.pprint(W)

detW, zeros_exact, poles_exact, zeros_numeric, poles_numeric = compute_det_and_zero_pole(W)
print("\n|W(s)| = det W(s) =")
sp.pprint(detW)

print("\nНули (точные, с кратностями):")
print(zeros_exact)
print("Полюса (точные, с кратностями):")
print(poles_exact)

print("\nНули (численно):")
print(zeros_numeric)
print("Полюса (численно):")
print(poles_numeric)

Ctr = controllability_matrix(A, B)
Obs = observability_matrix(A, C)
OutCtr = output_controllability_matrix(A, B, C)

rank_Ctr = int(Ctr.rank())
rank_Obs = int(Obs.rank())
rank_OutCtr = int(OutCtr.rank())

n = A.shape[0]
p = C.shape[0]

ctrl_full = rank_Ctr == n
obs_full = rank_Obs == n
outctrl_full = rank_OutCtr == p

stab_ok, bad_uncontrollable = is_stabilizable(A, B)
det_ok, bad_unobservable = is_detectable(A, C)

print("\n— Исследование свойств —")
print(f"Матрица управляемости:")
sp.pprint(Ctr)

print(f"Ранг матрицы управляемости = {rank_Ctr} из {n} → " + ("управляемая" if ctrl_full else "неуправляемая"))
print(f"Стабилизируемость: " + ("да" if stab_ok else f"нет; нестабильные неконтролируемые λ = {bad_uncontrollable}"))

print(f"Матрица наблюдаемости:")
sp.pprint(Obs)
print(f"Ранг матрицы наблюдаемости = {rank_Obs} из {n} → " + ("наблюдаемая" if obs_full else "ненаблюдаемая"))
print(f"Обнаруживаемость: " + ("да" if det_ok else f"нет; нестабильные ненаблюдаемые λ = {bad_unobservable}"))

print(f"Матрица управляемости по выходу:")
sp.pprint(OutCtr)
print(f"Ранг матрицы управляемости по выходу = {rank_OutCtr} из {p} → " + ("управляема по выходу" if outctrl_full else "НЕ управляема по выходу"))



# --- Временные характеристики ---
print("\n— Временные характеристики —")
try:
    t_sym, g_sym = symbolic_impulse_response(A, B, C)
    print("Весовая характеристика g(t) = C e^{At} B:")
    sp.pprint(g_sym)
except Exception as e:
    print("Не удалось получить g(t) символически:", e)

try:
    t_sym, h_sym = symbolic_step_response(A, B, C)
    print("Переходная характеристика h(t) для единичного входа:")
    sp.pprint(h_sym)
except Exception as e:
    print("Не удалось получить h(t) символически:", e)


# Численная визуализация g(t) и h(t)
t_grid = np.linspace(0, 10, 500)
try:
    g_vals = evaluate_matrix_function_t(g_sym, t_sym, t_grid)
    p, m, T = g_vals.shape
    for i in range(p):
        for j in range(m):
            plt.figure()
            plt.plot(t_grid, g_vals[i, j, :].real, label='Re')
            plt.plot(t_grid, g_vals[i, j, :].imag, '--', label='Im')
            plt.xlabel('t, s')
            plt.ylabel(f'g_{i+1}{j+1}(t)')
            plt.title('Весовая характеристика')
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.savefig(f'another_images/g_{i+1}{j+1}.png', dpi=150)
            plt.close()
except Exception as e:
    print("Не удалось построить g(t):", e)

try:
    h_vals = evaluate_matrix_function_t(h_sym, t_sym, t_grid)
    p, m, T = h_vals.shape
    for i in range(p):
        for j in range(m):
            plt.figure()
            plt.plot(t_grid, h_vals[i, j, :].real, label='Re')
            plt.plot(t_grid, h_vals[i, j, :].imag, '--', label='Im')
            plt.xlabel('t, s')
            plt.ylabel(f'h_{i+1}{j+1}(t)')
            plt.title('Переходная характеристика')
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.savefig(f'another_images/h_{i+1}{j+1}.png', dpi=150)
            plt.close()
except Exception as e:
    print("Не удалось построить h(t):", e)

# --- Частотные характеристики ---
print("\n— Частотные характеристики —")
w = np.linspace(0, 3, 1000)
Wjw = frequency_response(A, B, C, w)
p, m, _ = Wjw.shape

# АЧХ и ФЧХ
for i in range(p):
    for j in range(m):
        Wij = Wjw[i, j, :]
        mag = np.abs(Wij)
        phase = np.angle(Wij)
        # АЧХ/ФЧХ (ω в линейном масштабе)
        plt.figure()
        plt.subplot(2,1,1)
        plt.plot(w, mag)
        plt.grid(True)
        plt.ylabel('|W_{%d%d}(jω)|' % (i+1, j+1))
        plt.subplot(2,1,2)
        plt.plot(w, phase)
        plt.grid(True)
        plt.xlabel('ω, rad/s')
        plt.ylabel('arg W, rad')
        plt.tight_layout()
        plt.savefig(f'another_images/ACH_FCH_{i+1}{j+1}.png', dpi=150)
        plt.close()

        # ЛАЧХ/ЛФЧХ (Bode)
        plt.figure()
        plt.subplot(2,1,1)
        plt.semilogarithmic = True
        plt.semilogx(w, 20*np.log10(mag + 1e-12))
        plt.grid(True, which='both')
        plt.ylabel('20log10|W|, dB')
        plt.subplot(2,1,2)
        plt.semilogx(w, np.unwrap(phase))
        plt.grid(True, which='both')
        plt.xlabel('ω, rad/s')
        plt.ylabel('Phase, rad')
        plt.tight_layout()
        plt.savefig(f'another_images/Bode_{i+1}{j+1}.png', dpi=150)
        plt.close()
