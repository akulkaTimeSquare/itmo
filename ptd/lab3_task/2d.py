import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import sympy as sp

class RobotTracker:
    """
    Класс для моделирования слежения робота за траекторией на основе
    метода согласованного управления из предоставленного документа.
    """
    def __init__(self, params):
        self.params = params
        self.last_x_star = 0.0  # Начальное предположение для проекции

        # Инициализация символьных функций для траектории
        self._init_symbolic_functions()

    def _init_symbolic_functions(self):
        """
        Использует sympy для вычисления производных траектории и
        преобразования их в быстрые числовые функции.
        """
        x_sym = sp.Symbol('x')
        p = self.params
        A0, A1, wa, ws = p['A0'], p['A1'], p['wa'], p['ws']
        
        # Определение траектории g(x)
        g = (A0 + A1 * sp.cos(wa * x_sym)) * sp.sin(ws * x_sym)
        g_prime = sp.diff(g, x_sym)
        g_prime_prime = sp.diff(g_prime, x_sym)
        
        # Кривизна xi(x)
        xi = g_prime_prime / (1 + g_prime**2)**(3/2)
        
        # Производная кривизны по x
        dxi_dx = sp.diff(xi, x_sym)
        
        # Преобразование в числовые функции
        self.g_f = sp.lambdify(x_sym, g, 'numpy')
        self.g_prime_f = sp.lambdify(x_sym, g_prime, 'numpy')
        self.xi_f = sp.lambdify(x_sym, xi, 'numpy')
        self.dxi_dx_f = sp.lambdify(x_sym, dxi_dx, 'numpy')

    def _find_projection(self, x, y):
        """
        Находит точку (x_star, y_star) на траектории, ближайшую к (x, y).
        Это делается путем решения уравнения для перпендикуляра.
        """
        # Уравнение для поиска x_star: (x_star - x) + (g(x_star) - y) * g'(x_star) = 0
        func_to_solve = lambda xp: (xp - x) + (self.g_f(xp) - y) * self.g_prime_f(xp)
        
        # fsolve находит корень уравнения
        x_star = fsolve(func_to_solve, self.last_x_star)[0]
        self.last_x_star = x_star  # Сохраняем для следующего шага
        y_star = self.g_f(x_star)
        
        return x_star, y_star

    def dynamics(self, t, state):
        """
        Функция, описывающая динамику системы.
        Используется в решателе ОДУ (solve_ivp).
        """
        x, y, alpha, x_dot, y_dot, omega = state
        p = self.params

        # 1. Найти проекцию робота на траекторию
        x_star, y_star = self._find_projection(x, y)
        
        # 2. Вычислить геометрические свойства траектории в точке проекции
        g_prime = self.g_prime_f(x_star)
        alpha_star = np.arctan2(g_prime, 1)
        xi = self.xi_f(x_star)
        dxi_dx = self.dxi_dx_f(x_star)
        
        # 3. Рассчитать ошибки (e, delta) и их производные
        e = np.sign((y - y_star) - g_prime * (x - x_star)) * np.sqrt((x - x_star)**2 + (y - y_star)**2)
        
        delta = alpha - alpha_star
        delta = (delta + np.pi) % (2 * np.pi) - np.pi # Нормализация угла

        s_dot = x_dot * np.cos(alpha_star) + y_dot * np.sin(alpha_star)
        e_dot = -x_dot * np.sin(alpha_star) + y_dot * np.cos(alpha_star)
        
        delta_dot = omega - xi * s_dot
        xi_dot = dxi_dx * np.cos(alpha_star) * s_dot
        
        # 4. Вычислить управляющие воздействия (u_s, u_e, M)
        # Эти формулы основаны на принципе обратной линеаризации
        delta_v = s_dot - p['v_star']
        
        # Виртуальные управления (ускорения в системе координат Френе)
        u_s = xi * s_dot * e_dot - p['k_s'] * delta_v
        u_e = -xi * s_dot**2 - p['k_e1'] * e_dot - p['k_e2'] * e
        
        # Управляющий момент M
        s_ddot_approx = -p['k_s'] * delta_v
        M = p['J'] * (xi_dot * s_dot + xi * s_ddot_approx - p['k_d1'] * delta_dot - p['k_d2'] * delta)
        
        # 5. Преобразовать виртуальные управления в физические силы
        Fx_w = p['m'] * (u_s * np.cos(alpha_star) - u_e * np.sin(alpha_star))
        Fy_w = p['m'] * (u_s * np.sin(alpha_star) + u_e * np.cos(alpha_star))
        
        # 6. Вернуть производные состояния
        x_ddot = Fx_w / p['m']
        y_ddot = Fy_w / p['m']
        omega_dot = M / p['J']
        
        return [x_dot, y_dot, omega, x_ddot, y_ddot, omega_dot]

def run_simulation_and_plot():
    """
    Основная функция для запуска симуляции и построения графиков.
    """
    # Параметры из PDF и для симуляции
    params = {
        'm': 1.2, 'J': 1.0, 'v_star': 1.0,
        'x0': 0, 'y0': 1.5, 'alpha0': np.pi / 4,
        # Коэффициенты регуляторов (подобраны для устойчивости)
        'k_s': 10.0, 'k_e1': 20.0, 'k_e2': 20.0,
        'k_d1': 40.0, 'k_d2': 40.0,
        # Параметры траектории
        'A0': 1.0, 'A1': 1.0, 'wa': 2.0, 'ws': 1.0,
    }

    tracker = RobotTracker(params)

    # Начальные условия
    s0 = [params['x0'], params['y0'], params['alpha0'], 0, 0, 0]
    t_span = [0, 15]
    
    # Решение ОДУ
    sol = solve_ivp(tracker.dynamics, t_span, s0, dense_output=True, rtol=1e-6, atol=1e-6)
    
    # Обработка результатов
    t = np.linspace(t_span[0], t_span[1], 1000)
    S = sol.sol(t)
    x, y, alpha, _, _, _ = S

    # Пересчет ошибок для построения графика
    errors_e = []
    errors_delta = []
    errors_dv = []
    ss = []

    for i in range(len(t)):
        state_i = S[:, i]
        x_i, y_i, alpha_i, x_dot_i, y_dot_i, _ = state_i
        
        x_star, _ = tracker._find_projection(x_i, y_i)
        y_star = tracker.g_f(x_star)
        g_prime = tracker.g_prime_f(x_star)
        alpha_star = np.arctan2(g_prime, 1)

        e = np.sign((y_i - y_star) - g_prime * (x_i - x_star)) * np.sqrt((x_i - x_star)**2 + (y_i - y_star)**2)
        delta = alpha_i - alpha_star
        delta = (delta + np.pi) % (2 * np.pi) - np.pi
        
        s_dot = x_dot_i * np.cos(alpha_star) + y_dot_i * np.sin(alpha_star)
        delta_v = s_dot - tracker.params['v_star']

        ss.append(s_dot)
        errors_e.append(e)
        errors_delta.append(delta)
        errors_dv.append(delta_v)

    # Построение графиков
    
    # График 1: Ошибки
    plt.figure()
    plt.plot(t, errors_e, label='Ошибка отклонения')
    plt.plot(t, errors_delta, label='Ошибка поворота')
    plt.plot(t, errors_dv, label='Ошибка скорости')
    plt.title('Ошибки слежения за траекторией')
    plt.xlabel('t')
    plt.ylabel('f(t)')
    plt.legend()
    plt.grid(True)
    plt.savefig('errors.png')
    plt.show()

    # График 2: Траектории
    x_traj = np.linspace(min(x) - 1, max(x) + 1, 500)
    y_traj = tracker.g_f(x_traj)
    
    plt.figure()
    plt.plot(x_traj, y_traj, 'r--', label='Желаемая траектория')
    plt.plot(x, y, 'b', label='Фактическая траектория робота')
    plt.plot(params['x0'], params['y0'], 'go', markersize=10, label='Стартовая позиция')
    plt.title('Движение робота')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.grid(True)
    plt.savefig('track2d.png')
    plt.show()

if __name__ == '__main__':
    run_simulation_and_plot()