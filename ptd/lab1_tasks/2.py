import numpy as np
from scipy.integrate import quad

def rot_mat(phi):
    return np.array([[np.cos(phi), -np.sin(phi)],
                     [np.sin(phi),  np.cos(phi)]])

def parabola_length(x_end, k):
    # Интеграл длины дуги для y = k*x^3 => y' = 3*k*x^2
    f = lambda x: np.sqrt(1 + (3 * k * x**2)**2)
    length, _ = quad(f, 0, x_end)
    return length

points = np.array([
    [1, 1], [2, 1], [3, 2], [4, 2], [5, 3], 
    [6, 4], [7, 3], [8, 2], [9, 3], [9, 4], 
    [9, 5], [9, 6], [9, 7], [9, 8], [10, 9], [10, 10]
])

R = 0.8
k = 20 # Коэффициент кривизны параболы
total_length = 0
p_exit_last = points[0]

for i in range(1, len(points) - 1):
    p1, p2, p3 = points[i-1], points[i], points[i+1]
    
    # Расчет направлений и углов
    v1 = p2 - p1
    v2 = p3 - p2
    dir1 = v1 / np.linalg.norm(v1)
    dir2 = v2 / np.linalg.norm(v2)
    
    phi1 = np.arctan2(dir1[1], dir1[0])
    phi2 = np.arctan2(dir2[1], dir2[0])
    
    # Угол поворота (как в compute_angle)
    # В MATLAB det([dir1; dir2]) определяет сторону поворота
    sign_dir = np.sign(np.linalg.det(np.array([dir1, dir2])))
    
    # Внутренний угол sigma (corner_angle)
    dot_prod = np.clip(np.dot(-dir1, dir2), -1.0, 1.0)
    corner_angle = np.arccos(dot_prod)
    
    if abs(corner_angle - np.pi) < 1e-6:
        total_length += np.linalg.norm(p2 - p_exit_last)
        p_exit_last = p2
        continue

    d = R / np.abs(np.tan(corner_angle / 2))
    
    # Точки входа и выхода парабол
    parab1_entry = p2 - dir1 * d
    
    # x_end для параболы из формулы [1/(6*k*R); ...]
    xL_end = 1.0 / (6 * k * R)
    l_parab = parabola_length(xL_end, k)
    
    # Координаты конца первой параболы для вычисления h
    p1_exit_loc = np.array([xL_end, sign_dir * k * (xL_end**3)])
    parab1_exit = rot_mat(phi1) @ p1_exit_loc + parab1_entry
    
    parab2_exit = p2 + dir2 * d
    p2_entry_loc = np.array([xL_end, -sign_dir * k * (xL_end**3)])
    parab2_entry = rot_mat(phi2 - np.pi) @ p2_entry_loc + parab2_exit
    
    # Длина центральной дуги S3
    h = np.linalg.norm(parab2_entry - parab1_exit)
    # Угол дуги через хорду h
    arc_angle = 2 * np.arcsin(h / (2 * R))
    l_arc = R * arc_angle
    
    # Суммируем сегменты этого узла
    l_line = np.linalg.norm(parab1_entry - p_exit_last)
    
    total_length += + 2 * l_parab + l_arc
    p_exit_last = parab2_exit

# Последний прямой участок до финальной точки
total_length += np.linalg.norm(points[-1] - p_exit_last)

print(f"Итоговая длина: {total_length:.10f}")