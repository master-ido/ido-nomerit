import matplotlib.pyplot as plt
import numpy as np
from tools import Hermite, first_order_spline, cubic_spline_function, cubic_spline_derivative, cubic_spline_second_derivative

x_points_1 = [0, 2, 3, 4, 7, 8]
y_points_1 = [4, 2, 8, 10, 4, -2]
x_points_2 = [3, 2, 2.5, 4, 5, 4]
y_points_2 = [4, 3, 1, 2, 3.5, 4.5]
x_points_3 = [0, 1, 2, 3, 4, 5]
y_points_3 = [1, 1, 1, -1, -1, -1]

def step_function(x):
    if 0 <= x < 2.5:
        return 1
    if 2.5 <= x <=5:
        return -1
    return 0

def plot_spline_and_derivatives(x_points, y_points, num_points=100):
    spline = cubic_spline_function(x_points, y_points)
    spline_deriv = cubic_spline_derivative(x_points, y_points)
    spline_2nd_deriv = cubic_spline_second_derivative(x_points, y_points)
    plt.figure(figsize=(10, 6))
    for i in range(len(x_points) - 1):
        x_dense = np.linspace(x_points[i] + 1e-6, x_points[i + 1] - 1e-6, num_points)
        y_spline = [spline(x) for x in x_dense]
        y_deriv = [spline_deriv(x) for x in x_dense]
        y_2nd_deriv = [spline_2nd_deriv(x) for x in x_dense]
        plt.plot(x_dense, y_spline, color='blue', linewidth=4, label='Spline' if i == 0 else "")
        plt.plot(x_dense, y_deriv, color='green', linewidth=2, label="Spline'" if i == 0 else "")
        plt.plot(x_dense, y_2nd_deriv, color='red', linewidth=1, label="Spline''" if i == 0 else "")
    plt.scatter(x_points, y_points, color='black', label='Data points')
    plt.title('Cubic Spline and Its Derivatives')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_parametric_spline(x_points, y_points, num_points=100):
    x_points.extend([x_points[0], x_points[1]])
    y_points.extend([y_points[0], y_points[1]])
    t_points = [i for i in range(len(x_points))]
    spline_x = cubic_spline_function(t_points, x_points)
    spline_y = cubic_spline_function(t_points, y_points)
    plt.figure(figsize=(10, 6))
    for i in range(len(t_points) - 2):
        t_dense = np.linspace(t_points[i] + 1e-6, t_points[i + 1] - 1e-6, num_points)
        y_spline_x = [spline_x(t) for t in t_dense]
        y_spline_y = [spline_y(t) for t in t_dense]
        plt.plot(y_spline_x, y_spline_y, color='black', linestyle='-', label="" if i == 0 else "")
    plt.scatter(x_points, y_points, color='black', label='Data points')
    plt.title('Parametric Spline')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_cubic_first_order_hermite(x_points, y_points, num_points=100):
    cubic_spline = cubic_spline_function(x_points, y_points)
    first_order = first_order_spline(x_points, y_points)
    hermite_spline = Hermite(x_points, y_points, 0.5)
    plt.figure(figsize=(10, 6))
    for i in range(len(x_points) - 1):
        x_dense = np.linspace(x_points[i] + 1e-6, x_points[i + 1] - 1e-6, num_points)
        y_cubic_spline = [cubic_spline(x) for x in x_dense]
        y_first_order = [first_order(x) for x in x_dense]
        y_hermite_spline = [hermite_spline(x) for x in x_dense]
        y_step_function = [step_function(x) for x in x_dense]
        plt.plot(x_dense, y_cubic_spline, color='blue', label='Cubic' if i == 0 else "")
        plt.plot(x_dense, y_first_order, color='green', linewidth=5, label="First Order" if i == 0 else "")
        plt.plot(x_dense, y_hermite_spline, color='red', label="Hermite" if i == 0 else "")
        plt.plot(x_dense, y_step_function, c='black', label="Step Function" if i == 0 else "")
    plt.scatter(x_points, y_points, color='black', label='Data points')
    plt.title('Cubic Spline, First Order spline, Hermite')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    plot_spline_and_derivatives(x_points_1, y_points_1, 100)
    plot_parametric_spline(x_points_2, y_points_2, 100)
    plot_cubic_first_order_hermite(x_points_3, y_points_3, 100)