import numpy as np

A = np.array([[0, 1], [-1, 1]])
G = np.array([[-2, 2], [-2, -2]])
Y = np.array([[1, 1], [1, 1]])

V = np.vstack((Y, Y@G))

print("V = ", V)

B = np.array([[1, 2], [1, 0]])

print("BY = ", B@Y)

print(np.linalg.matrix_rank(B@Y))

b = np.array([[3], [1]])
h = np.array([[1, 1]])

print("U_A,b = ", np.hstack((b, A@b)))
print(np.linalg.matrix_rank(np.hstack((b, A@b))))

print("V_hG = ", np.vstack((h, h@G)))
print(np.linalg.matrix_rank(np.vstack((h, h@G))))
