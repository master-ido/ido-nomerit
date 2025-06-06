import numpy as np
from tools import finite_diff_order_4_error

def a(x):
    return (x ** 3) + (4 * x) - 15
ax, ah = 0, 0.25

def b(x):
    return (x ** 2) * np.cos(x)
bx, bh = 0.4, 0.1

def c(x):
    return np.tan(x / 3)
cx, ch = 3, 0.5

def d(x):
    return np.sin(0.5 * np.sqrt(x)) / x
dx, dh = 1, 0.2

def e(x):
    return np.e ** x + x
ex, eh = 2, 0.2

if __name__ == '__main__':
    forward_a, backward_a, centered_a = finite_diff_order_4_error(a, ax, ah)
    print(f'(a) Analytic result: 4, forward: {forward_a}, backward: {backward_a}, centered: {centered_a}')
    forward_b, backward_b, centered_b = finite_diff_order_4_error(b, bx, bh)
    print(f'(b) Analytic result: 0.6745418604, forward: {forward_b}, backward: {backward_b}, centered: {centered_b}')
    forward_c, backward_c, centered_c = finite_diff_order_4_error(c, cx, ch)
    print(f'(c) Analytic result: 1.141839607, forward: {forward_c}, backward: {backward_c}, centered: {centered_c}')
    forward_d, backward_d, centered_d = finite_diff_order_4_error(d, dx, dh)
    print(f'(d) Analytic result: -0.2600298981, forward: {forward_d}, backward: {backward_d}, centered: {centered_d}')
    forward_e, backward_e, centered_e = finite_diff_order_4_error(e, ex, eh)
    print(f'(e) Analytic result: 8.389056099, forward: {forward_e}, backward: {backward_e}, centered: {centered_e}')
