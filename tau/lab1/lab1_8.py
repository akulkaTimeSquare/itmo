import numpy as np

A = np.array([[7, -6, 9], [6, -5, 6], [-6, 3, -8]])
B = np.array([[1], [0], [0]])

#C = np.array([[0, -1, 1], [0, -3, 0]])
C = np.array([
    [6, 0, 6],
    [4, 1, 4]
])

D = np.array([[0], [0]])

U = np.hstack([B, A @ B, A @ A @ B])

Umat = np.hstack([C @ U, D])

print(Umat)
print(np.linalg.matrix_rank(Umat))