import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad

# Заданные матрицы
A = np.array([[7, -6,  9],
              [6, -5,  6],
              [-6, 3, -8]], dtype=float)
B = np.array([[2], [1], [-1]], dtype=float)

# Функция подынтегрального выражения
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

print(Wc)

pinvWc = np.linalg.pinv(Wc)
print(pinvWc@Wc@pinvWc)
