import numpy as np

# Определим матрицы A и B
A = np.array([
    [11, -2, 13],
    [6, -1, 6],
    [-6, -1, -8]
])
B = np.array([
    [2],
    [0],
    [0]
])

# 1. Проверка управляемости по Калману
def controllability_matrix(A, B):
    n = A.shape[0]
    ctrb = B
    for i in range(1, n):
        ctrb = np.hstack((ctrb, np.linalg.matrix_power(A, i) @ B))
    return ctrb

U = controllability_matrix(A, B)
rank_U = np.linalg.matrix_rank(U)
print("Матрица управляемости:\n", U)
print("Ранг матрицы управляемости:", rank_U)
if rank_U == A.shape[0]:
    print("Система полностью управляемая (по Калману)")
else:
    print("Система не полностью управляемая (по Калману)")

# 2. Найдём собственные числа матрицы A
eigvals = np.linalg.eigvals(A)
print("Собственные числа матрицы A:", eigvals)

# 3. Проверим каждое собственное число по критерию Хаутуса
for idx, lam in enumerate(eigvals):
    Hautus = np.hstack((A - lam * np.eye(A.shape[0]), B))
    rank_H = np.linalg.matrix_rank(Hautus)
    print(f"\nСобственное число λ{idx+1} = {lam:.4g}")
    print("Матрица Хаутуса:\n", Hautus)
    print("Ранг матрицы Хаутуса:", rank_H)
    if rank_H == A.shape[0]:
        print("Собственное число управляемо (по Хаутусу)")
    else:
        print("Собственное число не управляемо (по Хаутусу)")
