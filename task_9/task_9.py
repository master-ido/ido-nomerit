import math
import matplotlib.pyplot as plt
import numpy as np
from tools import present_finite_diff, present_shooting_with_RK4_2nd
def q_1(x):
    return 0.15
def p_1(x):
    return 0
def r_1(x):
    return 0
def dv_dt(x, v, y):
    return 0.15 * y
# present_finite_diff(p_1, q_1, r_1, 0, 240, 10, 150, 1)
present_shooting_with_RK4_2nd(dv_dt, 0, 240, 10, 150, 0.5)