from tools import plot_spline_and_derivatives, plot_parametric_spline, plot_cubic_first_order_hermite

x_points_1 = [0, 2, 3, 4, 7, 8]
y_points_1 = [4, 2, 8, 10, 4, -2]
x_points_2 = [3, 2, 2.5, 4, 5, 4]
y_points_2 = [4, 3, 1, 2, 3.5, 4.5]
x_points_3 = [0, 1, 2, 3, 4, 5]
y_points_3 = [1, 1, 1, -1, -1, -1]

if __name__ == '__main__':
    plot_spline_and_derivatives(x_points_1, y_points_1, 100)
    plot_parametric_spline(x_points_2, y_points_2, 100)
    plot_cubic_first_order_hermite(x_points_3, y_points_3, 100)