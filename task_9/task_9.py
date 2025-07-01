import math
import matplotlib.pyplot as plt
import numpy as np
from tools import present_finite_diff, present_shooting_with_RK4_2nd
def q_1(x):
    return 0.15
def q_2(x):
    return 0
def p_1(x):
    return 0
def p_2(x):
    return -12.5 / 70
def r_1(x):
    return 0
def r_2(x):
    return 9.81
def dv_dt_1(x, v, y):
    return 0.15 * y
def dv_dt_2(x, v, y):
    return 9.81 - (12.5 / 70) * v


present_finite_diff(p_1, q_1, r_1, 0, 240, 10, 150, 1)
present_shooting_with_RK4_2nd(dv_dt_1, 0, 240, 10, 150, 0.5)
present_finite_diff(p_2, q_2, r_2, 0, 500, 12, 0, 1)
present_shooting_with_RK4_2nd(dv_dt_2, 0, 500, 12, 0, 1)