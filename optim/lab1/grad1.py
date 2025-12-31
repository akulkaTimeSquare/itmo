import numpy as np
import matplotlib.pyplot as plt

# критерий и градиент
def J(xu):
    x, u = xu
    return 8*x**2 + 3*u**2 + 6*x*u + x + 4*u - 13

def grad_J(xu):
    x, u = xu
    return np.array([
        16*x + 6*u + 1,
        6*x + 6*u + 4
    ])

# параметры
gamma = 0.1
x = np.array([3.0, -4.0])

trajectory = [x]
J_values = [J(x)]

# итерации пока J не увеличивается
while True:
    x_new = x - gamma * grad_J(x)
    J_new = J(x_new)
    if J_new > J_values[-1]:
        break
    trajectory.append(x_new)
    J_values.append(J_new)
    x = x_new

trajectory = np.array(trajectory)
J_values = np.array(J_values)

# истинный минимум
x_star = 3/10
u_star = -29/30

# графики x(n), u(n)
n = np.arange(len(trajectory))

plt.figure()
plt.plot(n, trajectory[:,0], 'o-', label='x(n)')
plt.plot(n, trajectory[:,1], 's-', label='u(n)')

# горизонтальные линии на уровне минимума
plt.hlines(y=x_star, xmin=0, xmax=n[-1], linestyles='--', color='black')
plt.hlines(y=u_star, xmin=0, xmax=n[-1], linestyles='--', color='black')

# подписи линий сверху или снизу, чтобы не вылезали
plt.text(n[-1]-35, x_star + 0.25, r"$x^* = \frac{3}{10}$", fontsize=14, color='black', va='bottom')
plt.text(n[-1]-35, u_star - 0.25, r"$u^* = -\frac{29}{30}$", fontsize=14, color='black', va='top')

plt.xlabel('n')
plt.ylabel('f(n)')
plt.title(r"Состояния при методе наискорейшего спуска и $\gamma = 0.1$")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('images/grad_xu1.png')
plt.show()


# истинный минимум
xmin = np.array([3/10, -29/30])

# линии уровня
x_vals = np.linspace(-5, 5, 400)
u_vals = np.linspace(-5, 5, 400)
X, U = np.meshgrid(x_vals, u_vals)
Z = 8*X**2 + 3*U**2 + 6*X*U + X + 4*U - 13

plt.figure()
plt.contour(X, U, Z, levels=25)
plt.plot(trajectory[:,0], trajectory[:,1], 'o-')
plt.annotate(
    r"$\left(\frac{3}{10},-\frac{29}{30}\right)$",
    xy=(xmin[0], xmin[1]),
    xytext=(xmin[0]-1.5, xmin[1]+1),
    arrowprops=dict(arrowstyle="->"),
    fontsize=14,
)

plt.xlabel("x")
plt.ylabel("u")
plt.title(r"Траектория при методе наискорейшего спуска и $\gamma = 0.1$")
plt.grid(True)
plt.tight_layout()
plt.savefig('images/grad1.png')
plt.show()