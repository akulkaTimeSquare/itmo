import numpy as np
import matplotlib.pyplot as plt
import math

def func1(t, a, b):
    if abs(t) <= b:
        return a
    else:
        return 0

def furfunc1(w, a, b):
    return a*(2/math.pi)**(1/2)*np.sin(w*b)/w


a, b = 2, 6
y = np.vectorize(func1, otypes=[float])
t = np.linspace(-3*b, 3*b, 1000)
w = np.linspace(-3*b, 3*b, 1000)

ax1 = plt.subplot(211)
ax1.plot(t, y(t, a, b), label=f"Функция f(t)")
ax1.grid()
ax1.legend(loc="upper right")

ax2 = plt.subplot(212)
ax2.plot(w, furfunc1(w, a, b), label=f"Фурье-образ")
ax2.grid()
ax2.legend(loc="upper right")

plt.show()