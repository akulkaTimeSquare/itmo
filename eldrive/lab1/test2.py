import math

# Входные данные
d2 = 200            # Диаметр второго шкива, мм
omega0_ob_min = 1500  # Угловая скорость холостого хода, об/мин
omega_n_ob_min = 1440 # Номинальная угловая скорость, об/мин
Mn = 10             # Номинальный момент, Н·м
eta = 0.95          # Общее КПД системы
M2 = 21.507         # Момент на втором валу, Н·м
omega2 = 121.25     # Угловая скорость второго вала, рад/с

# Переводим об/мин в рад/с
omega0 = omega0_ob_min * 2 * math.pi / 60
omega_n = omega_n_ob_min * 2 * math.pi / 60

# Вычисляем жесткость характеристики h
h = Mn / (omega0 - omega_n)

# Вычисляем подкоренное выражение
numerator = 4 * M2 * omega2
denominator =  eta * h * (omega0 ** 2)
discriminant_inside = 1 - numerator / denominator

# Вычисляем корень
sqrt_term = math.sqrt(discriminant_inside)

# Вычисляем j
j_1 = (omega0 / (2 * omega2)) * (1 + sqrt_term)
j_2 = (omega0 / (2 * omega2)) * (1 - sqrt_term)

# Выбираем физически обоснованное значение j
j = j_1 if j_1 > 0 else j_2

# Вычисляем исправленный диаметр первого шкива
d1_tilde = d2 / j

# Вывод результатов с высокой точностью
print(f"omega0        = {omega0:.15f} рад/с")
print(f"omega_n       = {omega_n:.15f} рад/с")
print(f"h             = {h:.15f} Н·м·с/рад")
print(f"j             = {j:.15f}")
print(f"~d1           = {d1_tilde:.15f} мм")