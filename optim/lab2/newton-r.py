# Newton-Raphson method with annotated minimum
import numpy as np
import matplotlib.pyplot as plt

# Cost, gradient, Hessian
def J(xu):
    x, u = xu
    return 8*x**2 + 3*u**2 + 6*x*u + x + 4*u - 13

def grad_J(xu):
    x, u = xu
    return np.array([16*x + 6*u + 1,
                     6*x + 6*u + 4])

H = np.array([[16, 6],
              [6, 6]])
H_inv = np.linalg.inv(H)

# Newton iteration
x0 = np.array([3.0, -4.0])
trajectory = [x0]

x = x0.copy()
for _ in range(2):  # more steps are pointless for quadratic
    x = x - H_inv @ grad_J(x)
    trajectory.append(x)

trajectory = np.array(trajectory)
xmin = trajectory[-1]

# Contour plot
x_vals = np.linspace(-5, 5, 400)
u_vals = np.linspace(-5, 5, 400)
X, U = np.meshgrid(x_vals, u_vals)
Z = 8*X**2 + 3*U**2 + 6*X*U + X + 4*U - 13

plt.figure()
plt.contour(X, U, Z, levels=25)
plt.plot(trajectory[:, 0], trajectory[:, 1], 'o-', linewidth=2)

# Mark and annotate minimum
plt.scatter(xmin[0], xmin[1])
plt.annotate(
    r"$\left(\frac{3}{10},-\frac{29}{30}\right)$",
    xy=(xmin[0], xmin[1]),
    xytext=(xmin[0]-1.5, xmin[1]+1),
    arrowprops=dict(arrowstyle="->"),
    fontsize=14,
)

plt.xlabel("x")
plt.ylabel("u")
plt.title("Метод Ньютона-Рафсона")
plt.grid(True)
plt.tight_layout()
plt.savefig("images/nr.png")
plt.show()