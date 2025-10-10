import numpy as np

G1 = np.array([[-2, 1], [0, -2]])

Y1 = np.array([1, 1])

Q1 = np.array([Y1, Y1@G1])

print(np.linalg.matrix_rank(Q1))

G2 = np.array([[-20, 0], [0, -200]])

Y2 = np.array([1, 1])

Q2 = np.array([Y2, Y2@G2])

print(np.linalg.matrix_rank(Q2))