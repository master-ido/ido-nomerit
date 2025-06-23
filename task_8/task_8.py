import math
import matplotlib.pyplot as plt
import numpy as np
from tools import leapfrog, adams_4_with_rk_4, runga_kutta_4_order, runga_kutta_2_order, euler_method

def dm_dx(m,x):
    return 10 - 2 * x
def f(x, y):
    return (-2 * y) / (1 + x)
h = 0.5
y = 2
x = 0
final_x = 2

def d2x_dt2(y):
    return - (40 * y) / 2
initial_y = 0.7
initial_v = 0

# task 2:
def task_1(dy_dt, dv_dt, h):
