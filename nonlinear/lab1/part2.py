import os, math, numpy as np
import matplotlib.pyplot as plt

# Создаём директорию для сохранения изображений
outdir = 'images'
os.makedirs(outdir, exist_ok=True)

# Определяем системы
def sys1(t, x):
    x1,x2 = x
    return [-x1 + 2*x1**3 + x2, -x1 - x2]

def sys2(t, x):
    x1,x2 = x
    return [x1 + x1*x2, -x2 + x2**2 + x1*x2 - x1**3]

def sys3(t, x):
    x1,x2 = x
    return [x2, -x1 + x2*(1 - x1**2 + 0.1*x1**4)]

def sys4(t, x):
    x1,x2 = x
    r2 = x1**2 + x2**2
    return [(x1 - x2)*(1 - r2), (x1 + x2)*(1 - r2)]

def sys5(t, x):
    x1,x2 = x
    return [-x1**3 + x2, x1 - x2**3]

def sys6(t, x):
    x1,x2 = x
    return [-x1**3 + x2**3, x2**3 * x1 - x2**3]

systems = [
    (sys1, (-2.5,2.5,-2.5,2.5), 'system1'),
    (sys2, (-2.5,2.5,-2.5,2.5), 'system2'),
    (sys3, (-3,3,-3,3), 'system3'),
    (sys4, (-1.5,1.5,-1.5,1.5), 'system4'),  # ограничим область для корректного отображения
    (sys5, (-2,2,-2,2), 'system5'),
    (sys6, (-2,2,-2,2), 'system6')
]

# Равновесные точки
equilibria = {
    'system1': [(0,0), (1,-1), (-1,1)],
    'system2': [(0,0), (0,1), (1,-1)],
    'system3': [(0,0)],
    'system4': 'circle',  # окружность x1^2+x2^2=1 + центр
    'system5': [(0,0), (1,1), (-1,-1)],
    'system6': [(0,0), (1,1)]
}

# Построение фазовых портретов
for fun, bbox, name in systems:
    xmin,xmax,ymin,ymax = bbox
    x = np.linspace(xmin,xmax,25)
    y = np.linspace(ymin,ymax,25)
    X,Y = np.meshgrid(x,y)
    U = np.zeros_like(X)
    V = np.zeros_like(Y)

    # Заполнение векторного поля
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            u,v = fun(0, [X[i,j], Y[i,j]])
            if name=='system4' and X[i,j]**2 + Y[i,j]**2 > 1.5**2:
                U[i,j]=0
                V[i,j]=0
            else:
                U[i,j] = u
                V[i,j] = v

    fig, ax = plt.subplots(figsize=(6,6), dpi=150)
    ax.streamplot(X, Y, U, V, color='k', density=1.0)  # чёрные линии

    # Равновесные точки
    eq_pts = equilibria[name]
    if eq_pts == 'circle':
        theta = np.linspace(0,2*np.pi,200)
        ax.plot(np.cos(theta), np.sin(theta), 'r--')  # окружность пунктиром
        ax.plot(0,0,'ro')  # центр
    else:
        for ex,ey in eq_pts:
            ax.plot(ex, ey, 'ro')

    ax.set_xlim(xmin,xmax); ax.set_ylim(ymin,ymax)
    ax.set_xlabel('x1'); ax.set_ylabel('x2')
    ax.grid(True)
    plt.tight_layout()
    fname = os.path.join(outdir, f'{name}_eq.png')
    fig.savefig(fname)
    plt.close(fig)

print("Фазовые портреты сохранены в папке 'images/'.")
