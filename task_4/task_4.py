from tools import last, present_function, direct_method, lagrange, lagrange_polynom, present_polynom
import numpy as np

x = [-1, 0, 2]
y = [3, 1, 2]
x_1 = [-1, 0, 2, 3, 4]
y_1 = [0, 1, 9 ,25, 67]
x_2 = [-1, 0, 2]
y_2 = [0, 1, 9]
x_3 = [2, 3, 4]
y_3 = [9 ,25, 67]

f = lagrange_polynom(x_1, y_1)
g = lagrange_polynom(x_2, y_2)
h = lagrange_polynom(x_3, y_3)
def c(x):
    return np.cos(4 * np.pi * x)

if __name__ == '__main__':
    print(f"The result using Direct Method is {direct_method(x, y, -0.5)}", "\n")
    print(f"The result using Lagrange Method is {lagrange(x, y, -0.5)}", "\n")
    present_polynom(f, x_1, y_1)
    present_polynom(g, x_2, y_2)
    present_polynom(h, x_3, y_3)
    last(c)
