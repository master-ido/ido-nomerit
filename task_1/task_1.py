import numpy as np
from tools import present_function_3D, newton_raphson_2D, present_function, find_roots_with_synthetic_division, synthetic_division, find_roots_with_bisection, numeric_derivative, newton_raphson, find_roots_with_newton_raphson

def function(x):
    return (x ** 4) + (2 * x ** 3) - (7 * x ** 2) + 3

b_n = [1, 2, -7 , 0 ,3]

def f(x, y):
    return (4 * y ** 2) + (4 * y) - (52 * x) - 19

def g(x, y):
    return (169 * x ** 2) + (3 * y ** 2) - (111 * x) - (10 * y)

def function_3D(x, y):
    return np.sin(4 * y) * np.cos(0.5 * x)

if __name__ == '__main__':
    find_roots_with_bisection(function, -5, 3, 1)
    find_roots_with_newton_raphson(function, -5, 3, 1)
    find_roots_with_synthetic_division(b_n, -5, 3, 1)
    present_function(function)
    newton_raphson_2D(f, g, -0.01, -0.01)
    present_function_3D(function_3D)







