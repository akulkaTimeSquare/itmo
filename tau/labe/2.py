import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

def get_eigenvalues(A):
    return np.linalg.eig(A)[0]

A = np.array([[0, 1], [-1, 1]])
print(get_eigenvalues(A))

B = np.array([[1, 2], [1, 0]])
C = np.array([[1, -2], [0, 3]])
D = np.array([[-3, 0], [0, 1]])
Cz = np.array([[1, 2], [4, 0]])
Dz = np.array([[-3, 0], [0, 1]])

Bf = np.array([[1, 2], [1, 3]])
Df = np.array([[1, 0], [0, 1]])

Vy = np.array(np.vstack((C, C@A)))
print("Vy:", Vy)
print("rank(Vy):", np.linalg.matrix_rank(Vy))

Vz = np.array(np.vstack((Cz, Cz@A)))
print("\nVz:", Vz)
print("rank(Vz):", np.linalg.matrix_rank(Vz))

Uy = np.array(np.hstack((C@B, C@A@B, D)))
print("\nUy:", Uy)
print("rank(Uy):", np.linalg.matrix_rank(Uy))

Uz = np.array(np.hstack((Cz@B, Cz@A@B, Dz)))
print("\nUz:", Uz)
print("rank(Uz):", np.linalg.matrix_rank(Uz))


def compute_transfer_matrix(A, B, C, D):
    s = sp.symbols('s')
    A_sp = sp.Matrix(A.tolist())
    B_sp = sp.Matrix(B.tolist())
    C_sp = sp.Matrix(C.tolist())
    D_sp = sp.Matrix(D.tolist())
    I = sp.eye(A_sp.shape[0])
    W = C_sp * (s*I - A_sp).inv() * B_sp + D_sp
    return sp.simplify(W)

Wy = compute_transfer_matrix(A, B, C, D)
print("\nWy:")
sp.pprint(Wy)
Wz = compute_transfer_matrix(A, B, Cz, Dz)
print("\nWz:")
sp.pprint(Wz)


# Determinants of Wy and Wz
det_Wy = sp.simplify(sp.factor(Wy.det()))
det_Wz = sp.simplify(sp.factor(Wz.det()))

print("\n det(Wy) =")
sp.pprint(det_Wy)

print("\n det(Wz) =")
sp.pprint(det_Wz)

