import math
import matplotlib.pyplot as plt
import numpy as np
from tools import present_function, present_leapfrog, present_adams_4_with_rk_4, present_rk_4, present_rk_2, present_euler, present_RK4_2nd, runga_kutta_4_order_2nd, leapfrog, adams_4_with_rk_4, runga_kutta_4th_order, runga_kutta_2nd_order, euler_method, present_euler_2nd

def dv_dt(t, v, y):
    return -20 * y
def analytic_func(t):
    return 0.7 * np.cos((np.sqrt(20) * t))

def dm_dx(x, m):
    return 10 - 2 * x
def f(x, y):
    return (-2 * y) / (1 + x)
def m(x):
    return 10 * x - x ** 2
if __name__ == '__main__':
    # 1. a) euler with h = 0.1:
    present_euler_2nd(dv_dt, 0, 0, 0.7, 2.5, 0.1)
    # 1. a) RK4 with h = 0.1:
    present_RK4_2nd(dv_dt, 0, 0, 0.7, 2.5, 0.1)
    # 1. b) euler with h = 0.05:
    present_euler_2nd(dv_dt, 0, 0, 0.7, 2.5, 0.05)
    # 1. c) leapfrog with h = 0.1
    present_leapfrog(dv_dt, 0, 0, 0.7, 0.1, 2.5)
    # analytic solution:
    present_function(analytic_func, 0, 2.5)
    # 2. a) Euler with h = 0.25
    present_euler(dm_dx, 0, 0, 10, 0.25)
    # 2. a) Euler with h = 0.05
    present_euler(dm_dx, 0, 0, 10, 0.05)
    # 2. b) RK2 with h = 0.25
    present_rk_2(dm_dx, 0, 10, 10, 0.25)
    # 2. c) analytic solution:
    present_function(m, 0, 10)
    # 3. Adam's method with RK4:
    present_adams_4_with_rk_4(f, 0, 2, 10, 0.5)



