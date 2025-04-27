from tools import present_function, direct_method, lagrange, lagrange_polynom, present_polynom
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

def last(function):
    chosen_number = float(input("Choose a number between 0 to 0.5: "))
    print("\n" f"The analytic result is {round(function(chosen_number), 7)}" "\n")
    measure_points_1 = np.linspace(0, 0.5, 4)
    measure_points_2 = np.linspace(0, 0.5, 8)
    measure_results_1 = function(measure_points_1)
    measure_results_2 = function(measure_points_2)
    print(f"The interpolation result using 4 measure points is {round(lagrange(measure_points_1, measure_results_1, chosen_number), 7)}" "\n")
    present_polynom(function, measure_points_1, measure_results_1)
    print(f"The interpolation result using 8 measure points is {round(lagrange(measure_points_2, measure_results_2, chosen_number), 7)}")
    present_polynom(function, measure_points_2, measure_results_2)
    present_function(function)

if __name__ == '__main__':
    print(f"The result using Direct Method is {round(direct_method(x, y, -0.5), 7)}", "\n")
    print(f"The result using Lagrange Method is {round(lagrange(x, y, -0.5), 7)}", "\n")
    present_polynom(f, x_1, y_1)
    present_polynom(g, x_2, y_2)
    present_polynom(h, x_3, y_3)
    last(c)
