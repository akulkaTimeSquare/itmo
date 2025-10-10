import numpy as np
import matplotlib.pyplot as plt

def x(t):
    if -2 <= t < -1:
        return 0.5
    elif -1 <= t < 0:
        return -t-0.5
    elif 0 <= t < 1.5:
        return -t-0.5
    elif 1.5 <= t < 2:
        return t-3.5
    elif 2 <= t < 4:
        return t-3.5
    elif 4 <= t < 5.5:
        return t-3.5
    elif 5.5 <= t <= 6:
        return -t+7.5

ran = np.arange(-2, 6, 0.01)
y = []
for t in ran:
    y += [x(t)]
plt.plot(ran, y, label="$x(t)$")
plt.grid(True)
plt.legend()
plt.show()