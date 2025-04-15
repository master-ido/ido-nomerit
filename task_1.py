import numpy as np

def bisection(function, a, b):
    if abs(a - b) < 10 ** -6:
        return round((a + b) / 2, 7)
    u = function(a)
    v = function(b)
    if u * v == 0:
        if u == 0:
            return a
        else:
            return b
    if u * v < 0:
        c = (a + b) / 2
        w = function(c)
        if u * w < 0:
            return bisection(function, a, c)
        else:
            return bisection(function, c, b)

def function(x):
    return (x ** 4) + (2 * x ** 3) - (7 * x ** 2) + 3

def find_roots_with_bisection(function, start, end, step):
    # use different steps, to make sure you didnt miss a root
    for i in np.arange(start, end, step):
        if bisection(function, i, i + step) == None:
            print(f"There are no roots between {i}, {i + step}")
        else:
            print(f"The root between {i}, {i + step} is {bisection(function, i, i + step)}")

find_roots_with_bisection(function, -5, 3, 1)


