import numpy as np
import matplotlib.pyplot as plt

dt = 0.05
t = np.arange(-100, 100, dt)
y = np.sin(t)
a = 0.1
y = y + a*(np.random.rand(len(t))-0.5)
num = [0]
for k in range(1, len(y)):
    num += [(y[k] - y[k-1])/dt]
num = np.array(num)

furie = []
vx = np.arange(-100, 100, dt)
for v in vx:
    func = y*np.exp(-2*np.pi*1j*v*t)
    im = np.trapz(func, t, dt)
    furie += [im]
furie = np.array(furie)
der_furie = 2*np.pi*1j*vx*np.array(furie)

der = []
for k in t:
    der += [np.trapz(der_furie*np.exp(2*np.pi*1j*vx*k), vx)]

plt.plot(vx, np.real(der))
"""
# plt.plot(t, num)
plt.plot(t, np.cos(t))
plt.plot(t, der)
plt.ylim(-2, 2)
plt.grid(True)
"""
"""
#plt.plot(t, num)
fig, axs = plt.subplots(2)
axs[0].plot(t, np.real(der))
axs[0].plot(t, num)
# axs[0].plot(t, np.cos(t))
axs[0].set_ylim(-2, 2)
axs[0].grid(True)
axs[1].plot(t, np.imag(der))
axs[1].grid(True)
"""
plt.show()