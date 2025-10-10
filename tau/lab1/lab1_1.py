import numpy as np

A = np.array([[7, -6, 9], [6, -5, 6], [-6, 3, -8]])
B = np.array([[2], [1], [-1]])

x1d = np.array([[-2], [1], [-1]])
x1dd = np.array([[-5], [4], [-1]])

U = np.hstack([B, A @ B, A @ A @ B])


Ud = np.hstack([B, A @ B, A @ A @ B, x1d])

Udd = np.hstack([B, A @ B, A @ A @ B, x1dd])

print(np.linalg.matrix_rank(U))
print(np.linalg.matrix_rank(Ud))
print(np.linalg.matrix_rank(Udd))

x1 = x1dd

eigvals = np.linalg.eigvals(A)

#print("Eigenvalues:", eigvals)

H0 = np.hstack([A - eigvals[0] * np.eye(3), B])
H1 = np.hstack([A - eigvals[1] * np.eye(3), B])
H2 = np.hstack([A - eigvals[2] * np.eye(3), B])

#print(H2)

#print(np.linalg.matrix_rank(H0))
#print(np.linalg.matrix_rank(H1))
#print(np.linalg.matrix_rank(H2))

P = np.array([[-1, -1.5, -0.5], [0, -1, 0], [1, 1, 0]])
invP = np.linalg.inv(P)
D = np.array([[-2, 0, 0], [0, -2, 3], [0, -3, -2]])
#print(invP @ B)