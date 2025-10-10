import numpy as np
import matplotlib.pyplot as plt

def f(t):
    t -= t//T*T
    if t >= t0 and t < t1: return a
    return b

a = 1
b = 2
t0 = 1
t1 = 2
t2 = 3    
dt = 1e-5
T = t2-t0
t = np.arange(-5, 5, dt)
#t = np.arange(-3*T, 3*T, dt)
func = []
for i in t:
    func += [f(i)]
func = np.array(func)

N = 2
coef_t = np.arange(t0, t2, dt)
y = []
for i in coef_t:
    y += [f(i)]
y = np.array(y)
an = [np.trapz(2/T*y, coef_t)]
cn = [np.trapz(1/T*y, coef_t)]
bn = [0]
g = an[0]/2
cg = cn[0]
coef_g = an[0]/2
coef_cg = cn[0]
for n in range(1, N+1):
    an += [np.trapz(2/T*y*np.cos(2*np.pi*n/T*coef_t), coef_t)]
    bn += [np.trapz(2/T*y*np.sin(2*np.pi*n/T*coef_t), coef_t)]
    cn = [np.trapz(1/T*y*np.exp(1j*2*np.pi*n/T*coef_t), coef_t)] +\
        cn + [np.trapz(1/T*y*np.exp(-1j*2*np.pi*n/T*coef_t), coef_t)]
    g += an[-1]*np.cos(2*np.pi*n/T*t) + bn[-1]*np.sin(2*np.pi*n/T*t)
    coef_g += an[-1]*np.cos(np.pi*n*coef_t) + bn[-1]*np.sin(np.pi*n*coef_t)
    cg += cn[0]*np.exp(-1j*2*np.pi*n/T*t) + cn[-1]*np.exp(1j*2*np.pi*n/T*t)
    coef_cg += cn[0]*np.exp(-1j*2*np.pi*n/T*coef_t) + cn[-1]*np.exp(1j*2*np.pi*n/T*coef_t)

print(f"an: {list(map(lambda x: np.round(x, 3), an))}")
print(f"bn: {list(map(lambda x: np.round(x, 3), bn))}")
print(f"cn: {list(map(lambda x: np.round(x, 3), cn))}")

#plt.plot(t, func)
#plt.plot(t, g)

#plt.plot(coef_t, y, label="Квадратная волна")
#plt.plot(coef_t, coef_g, label="Вещественный ряд")

plt.plot(coef_t, y, label="Квадратная волна")
plt.plot(coef_t, coef_cg, label="Комплексный ряд")

plt.grid(True)
#plt.ylim(top=2.15)
#plt.legend()
#plt.savefig(f"Documents/Frequency methods/Lab1/images/fig.png", bbox_inches='tight')
#plt.savefig(f"Documents/Frequency methods/Lab1/images/k={N}.png", bbox_inches='tight')
plt.show()