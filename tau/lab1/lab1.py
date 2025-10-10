import numpy as np

A = np.array([[7, -6, 9], [6, -5, 6], [-6, 3, -8]])
B = np.array([[-2], [-1], [2]])
x1 = np.array([[-5], [-3], [3]])

U = np.hstack([B, A @ B, A @ A @ B])

eigvals = np.linalg.eigvals(A)

H0 = np.hstack([A - eigvals[0] * np.eye(3), B])
H1 = np.hstack([A - eigvals[1] * np.eye(3), B])
H2 = np.hstack([A - eigvals[1] * np.eye(3), B])


P = np.array([[-1, -1.5, -0.5], [0, -1, 0], [1, 1, 0]])
invP = np.linalg.inv(P)
D = np.array([[-2, 0, 0], [0, -2, 3], [0, -3, -2]])

print()
