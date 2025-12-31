import numpy as np

points = np.array([
    [1, 1], [2, 1], [3, 2], [4, 2], [5, 3], 
    [6, 4], [7, 3], [8, 2], [9, 3], [9, 4], 
    [9, 5], [9, 6], [9, 7], [9, 8], [10, 9], [10, 10]
])

R = 0.8
total_length = 0

def length(v):
    return np.sqrt(v[0]**2 + v[1]**2)

# Проходим по всем сегментам
# Логика: прибавляем полную длину сегмента, затем корректируем углы
segments = []
for i in range(len(points) - 1):
    dist = length(points[i+1] - points[i])
    total_length += dist

# Корректировка на углах (вычитаем 2*d, прибавляем дугу)
for i in range(1, len(points) - 1):
    p_prev = points[i-1]
    p_curr = points[i]
    p_next = points[i+1]
    
    v1 = p_curr - p_prev
    v2 = p_next - p_curr
    
    # Угол между векторами
    angle = np.arctan2(v2[1], v2[0]) - np.arctan2(v1[1], v1[0])
    
    # Нормализация угла (-pi, pi)
    angle = (angle + np.pi) % (2 * np.pi) - np.pi
    
    if abs(angle) > 1e-6: # Если есть поворот
        theta = abs(angle)
        d = R * np.tan(theta / 2)
        arc = R * theta
        
        # Коррекция: убираем прямые куски (d с каждой стороны), добавляем дугу
        total_length = total_length - 2*d + arc

print(f"Итоговая длина: {total_length:.10f}")