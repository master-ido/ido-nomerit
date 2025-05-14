import math
import matplotlib.pyplot as plt
import numpy as np

# Reading a matrix
def txt_to_matrix(txt):
    with open(txt, "r") as file:
        lines = file.readlines()
        matrix = []
        for line in lines:
            str_row = line.split()
            row = list(map(float, str_row))
            matrix.append(row)
    return matrix

# Reading a vector
def txt_to_vector(txt):
    with open("vec_input.txt", "r") as v_file:
        vec = v_file.readlines()
        vector = [float(vec[i]) for i in range(len(vec))]
    return vector

# presentable matrix
def matrix_print(matrix):
    for row in matrix:
        print(", ".join(map(str, row)))

# minor for the determinant
def get_minor(matrix, row, col):
    minor = []
    for i in range(len(matrix)):
        if i == row:
            continue
        minor_row = []
        for j in range(len(matrix[i])):
            if j == col:
                continue
            minor_row.append(matrix[i][j])
        minor.append(minor_row)
    return minor

# determinant
def determinant(matrix):
    if len(matrix) == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    det = 0
    for col in range(len(matrix)):
        minor = get_minor(matrix, 0, col)
        det += ((-1) ** col) * matrix[0][col] * determinant(minor)
    return det

# trace
def get_trace(matrix):
    trace = 0
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if i == j:
                trace += matrix[i][j]
    return(trace)

# multiplying matrix with a vector
def multiplying(matrix, vector):
    result = []
    for i in range(len(matrix)):
        vec_line = 0
        for j in range(len(matrix[i])):
            vec_line += matrix[i][j] * vector[j]
        result.append(vec_line)
    return result

# NMO
def my_nmo_base(x, t_0, v):
    t_squared = t_0 + (x ** 2)/(v ** 2)
    return math.sqrt(t_squared)

def vector_to_txt(vector):
    with open("mvm_out.txt", "w") as r_file:
        for i in vector:
            r_file.write(str(i) + "\n")

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

def find_roots_with_bisection(function, start, end, step):
    # use different steps, to make sure you didnt miss a root
    roots = []
    for i in np.arange(start, end, step):
        if bisection(function, i, i + step) == None:
            continue
            # print(f"There are no roots between {i}, {i + step}")
        if bisection(function, i, i + step) in roots:
            continue
        roots.append(round(bisection(function, i, i + step), 5))
    print(f"Using Bisection, the roots between {start}, {end} are {roots}")

def numeric_derivative(function, x):
    return (function(x + 10 ** -7) - function(x)) / (10 ** -7)


def newton_raphson(function, x):
    x_new = x - (function(x) / numeric_derivative(function, x))
    if abs(x_new - x) < 10 ** -7:
        return x_new
    return newton_raphson(function, x_new)

def find_roots_with_newton_raphson(function, start, end, step):
    roots = []
    for i in np.arange(start, end, step):
        if round(newton_raphson(function, i),5) in roots:
            continue
        roots.append(round(newton_raphson(function, i), 5))
    print(f'Using Newton-Raphson, the roots between {start}, {end} are {roots}')


def synthetic_division(b_n, x):
    # creating lists
    c_n = []
    d_n = []
    c_n.append(b_n[0])
    for i in range(1, len(b_n) - 1):
        c_n.append(b_n[i] + x * c_n[i - 1])
    # c_n is ready, thus r_0 is:
    r_0 = b_n[-1] + x * c_n[-1]
    d_n.append(c_n[0])
    for i in range(1, len(c_n) - 1):
        d_n.append(c_n[i] + x * d_n[i - 1])
    # d_n is ready, thus r_1 is:
    r_1 = c_n[-1] + x * d_n[-1]
    # newton raphson with our r's:
    x_new = x - (r_0 / r_1)
    if abs(x_new - x) < 10 ** -7:
        return x_new
    else:
        return synthetic_division(b_n, x_new)

def find_roots_with_synthetic_division(b_n, start, end, step):
    roots = []
    for i in range(start, end, step):
        if i == 0:
            continue
        if round(synthetic_division(b_n, i), 5) in roots:
            continue
        roots.append(round(synthetic_division(b_n, i), 5))
    print(f'Using Synthetic Division, the roots between {start}, {end} are {roots}')


def present_function(function):
    x = np.linspace(0, 0.5, 500)
    y = function(x)
    plt.plot(x, y, color='blue')
    plt.axhline(0, color='gray', linewidth=0.5)  # x-axis
    plt.axvline(0, color='gray', linewidth=0.5)  # y-axis
    plt.title("Plot of f(x)")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    # plt.legend()
    plt.grid(True)
    plt.show()


def partial_derivative_dx(function, x, y):
    df_dx = (function(x + 1e-5, y) - function(x, y)) / 1e-5
    return df_dx

def partial_derivative_dy(function, x, y):
    df_dy = (function(x, y + 1e-5) - function(x, y)) / 1e-5
    return df_dy

def newton_raphson_2D(f, g, x, y, x0=None, y0=None):
    if x0 is None and y0 is None:
        x0, y0 = x, y

    df_dx = partial_derivative_dx(f, x, y)
    df_dy = partial_derivative_dy(f, x, y)
    dg_dx = partial_derivative_dx(g, x, y)
    dg_dy = partial_derivative_dy(g, x, y)

    x_new = x - ((f(x, y) * dg_dy) - g(x, y) * df_dy) / (df_dx * dg_dy - dg_dx * df_dy)
    y_new = y - ((g(x, y) * df_dx) - f(x, y) * dg_dx) / (df_dx * dg_dy - dg_dx * df_dy)

    if abs(x_new - x) and abs(y_new - y) < 10 ** -7:
        print(f'Using Newton-Raphson, the solution for starting point ({x0}, {y0}) is ({round(x_new, 5)}, {round(y_new, 5)})')
        return
    return newton_raphson_2D(f, g, x_new, y_new, x0, y0)



def present_function_3D(function):
    x = np.linspace(-10, 10, 50)
    y = np.linspace(-5, 5, 50)
    X, Y = np.meshgrid(x, y)
    Z = function(X, Y)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='Purples', edgecolor='none')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.show()

def partial_pivoting(matrix, vector):
    for j in range(len(matrix) - 1):
        max_row = j
        max_val = abs(matrix[j][j])
        for i in range(j + 1, len(matrix)):
            if abs(matrix[i][j]) > max_val:
                max_val = abs(matrix[i][j])
                max_row = i
        if max_row != j:
            matrix[j], matrix[max_row] = matrix[max_row], matrix[j]
            vector[j], vector[max_row] = vector[max_row], vector[j]
    return matrix, vector

def gauss(matrix, vector):
    matrix, vector = partial_pivoting(matrix, vector)
    for i in range(len(matrix)):
        pivot = matrix[i][i]
        for j in range(i, len(matrix)):
            matrix[i][j] = matrix[i][j] / pivot
        vector[i] = vector[i] / pivot
        for k in range(i + 1, len(matrix)):
            factor = matrix[k][i]
            for j in range(i, len(matrix)):
                matrix[k][j] = matrix[k][j] - (factor * matrix[i][j])
            vector[k] = vector[k] - (factor * vector[i])
    result = [0 for _ in range(len(matrix))]
    for i in range(len(matrix) - 1, -1, -1):
        result[i] = vector[i] - sum(matrix[i][j] * result[j] for j in range(i + 1, len(matrix)))
    return [round(x, 5) for x in result]

def gauss_seidel(matrix, vector, x=None):
    if x is None:
        x = np.random.rand(len(matrix))
    matrix, vector = partial_pivoting(matrix, vector)
    reference_x = x.copy()
    for i in range(len(matrix)):
        sum = 0
        for j in range(len(matrix)):
            if j == i:
                continue
            sum += matrix[i][j] * x[j]
        x[i] = (vector[i] - sum) / matrix[i][i]
    for i in range(len(x)):
        if abs(x[i] - reference_x[i]) > (10 ** -5):
            return gauss_seidel(matrix, vector, x)
        return [round(i, 3) for i in x]

def LU_Decomposition(matrix):
    L = np.identity(len(matrix))
    U = np.zeros([len(matrix), len(matrix)])
    for j in range(len(matrix)):
        for i in range(len(matrix)):
            if i <= j:
                U[i][j] = matrix[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
            if i > j:
                L[i][j] = (matrix[i][j] - sum(L[i][k] * U[k][j] for k in range(j))) / U[j][j]
    return L, U

def solve_with_LU(matrix, vector):
    L, U = LU_Decomposition(matrix)
    y = [vector[0] / L[0][0] for _ in range(len(matrix))]
    for i in range(1, len(matrix)):
        y[i] = (vector[i] - sum(L[i][j] * y[j] for j in range(i))) / L[i][i]
    result = [y[-1] / U[-1][-1] for _ in range(len(matrix))]
    for i in range(len(matrix) - 1, -1, -1):
        result[i] = (y[i] - sum(U[i][j] * result[j] for j in range(i + 1, len(matrix)))) / U[i][i]
    return [round(x, 4) for x in result]

def inverse_matrix(matrix):
    inv_mat = np.zeros([len(matrix), len(matrix)])
    for i in range(len(matrix)):
        I = [0 for _ in range(len(matrix))]
        I[i] = 1
        x = solve_with_LU(matrix, I)
        for k in range(len(matrix)):
            inv_mat[k][i] = x[k]
    return inv_mat

def multiplying_matrix(matrix_1, matrix_2):
    result = np.zeros([len(matrix_1), len(matrix_1)])
    for i in range(len(matrix_1)):
        for j in range(len(matrix_1)):
            for k in range(len(matrix_1)):
                result[i][j] += matrix_1[i][k] * matrix_2[k][j]
    return result

def direct_method(x_points, y_points, x):
    matrix = np.ones((len(x_points), len(x_points)))
    for i in range(len(x_points)):
        for j in range(1, len(x_points)):
            matrix[i][j] = x_points[i] ** j
    a = solve_with_LU(matrix, y_points)
    result = a[0]
    for i in range(1, len(a)):
        result += a[i] * (x ** i)
    return result

def lagrange_l_polynom(x_points, i, x):
    l = 1
    for j in range(len(x_points)):
        if j == i:
            continue
        l *= (x - x_points[j]) / (x_points[i] - x_points[j])
    return l


def lagrange(x_points, y_points, x):
    result = 0
    for i in range(len(x_points)):
        result += lagrange_l_polynom(x_points, i, x) * y_points[i]
    return result


def lagrange_polynom(x_points, y_points):
    return lambda x: lagrange(x_points, y_points, x)


def present_polynom(function, x_points, y_points):
    x = np.linspace(min(x_points) - 0.01, max(x_points) + 0.01, 500)
    f_x = [function(p) for p in x]
    plt.plot(x, y_points, color='blue', label='Lagrange interpolation')
    plt.scatter(x_points, y_points, color='red', label='Original Points')
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.axvline(0, color='gray', linewidth=0.5)
    plt.title("Lagrange Interpolation Polynomial")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.legend()
    plt.grid(True)
    plt.show()

def banded_system(mat, res):
    for i in range(1, len(mat)):
        mat[i][i] -= mat[i - 1][i] * (mat[i][i - 1] / mat[i - 1][i - 1])
        res[i] -= res[i - 1] * (mat[i][i - 1] / mat[i - 1][i - 1])
        mat[i][i - 1] = 0
    x = [res[i] / mat[i][i] for i in range(len(mat))]
    for i in range(len(mat) - 2, -1, -1):
        x[i] = (res[i] - mat[i][i + 1] * x[i + 1]) / mat[i][i]
    return [round(i, 7) for i in x]

def cubic_spline_parameters(x_points, y_points):
    h = [x_points[i + 1] - x_points[i] for i in range(len(x_points) - 1)]
    u = [2 * (h[i] + h[i - 1]) for i in range(1, len(h))] # the length might be too long
    k = [(y_points[i + 1] - y_points[i]) / h[i] for i in range(len(h))]
    v = [6 * (k[i] - k[i - 1]) for i in range(1, len(k))]
    return h, u, v


def m_cubic_splines(x_points, y_points):
    h, u, v = cubic_spline_parameters(x_points, y_points)
    mat = np.zeros((len(u), len(u)))
    for i in range(len(u)):
        mat[i][i] = u[i]
    for i in range(len(u) - 1):
        mat[i][i + 1] = h[i + 1]
        mat[i + 1][i] = h[i + 1]
    return [0] + banded_system(mat, v) + [0], h
def cubic_spline_coeff(x_points, y_points):
    m, h = m_cubic_splines(x_points, y_points)
    a = [(m[i + 1] - m[i]) / (6 * h[i]) for i in range(len(h))]
    b = [(m[i] * x_points[i + 1] - m[i + 1] * x_points[i]) / (2 * h[i]) for i in range(len(m) - 1)]
    c = [((m[i + 1] * 3 * (x_points[i] ** 2) - m[i] * 3 * (x_points[i + 1] ** 2)) + (6 * (y_points[i + 1] - y_points[i])) + ((h[i] ** 2) * (m[i] - m[i + 1]))) / (6 * h[i]) for i in range(len(m) - 1)]
    d = [(m[i] * (x_points[i + 1] ** 3) - m[i + 1] * (x_points[i] ** 3) + 6 * (x_points[i + 1] * y_points[i] - x_points[i] * y_points[i + 1]) + (h[i] ** 2) * (m[i + 1] * x_points[i] - m[i] * x_points[i + 1])) / (6 * h[i]) for i in range(len(m) - 1)]
    s = []
    for i in range(len(a)):
        s_use = [a[i], b[i], c[i], d[i]]
        s.append(s_use)
    return s

def cubic_spline_function_result(x_points, y_points, x):
    s = cubic_spline_coeff(x_points, y_points)
    for i in range(len(x_points) - 1):
        if x_points[i] < x < x_points[i + 1]:
            result = 0
            for j in range(len(s[0])):
                result += s[i][j] * (x ** (len(s[0]) - 1 - j))
            return result

def cubic_spline_function(x_points, y_points):
    return lambda x: cubic_spline_function_result(x_points, y_points, x)

def cubic_spline_derivative_result(x_points, y_points, x):
    s = cubic_spline_coeff(x_points, y_points)
    for i in range(len(x_points) - 1):
        if x_points[i] < x < x_points[i + 1]:
            result = (3 * s[i][0] * (x ** 2)) + (2 * s[i][1] * x) + s[i][2]
            return result

def cubic_spline_derivative(x_points, y_points):
    return lambda x: cubic_spline_derivative_result(x_points, y_points, x)

def cubic_spline_second_derivative_result(x_points, y_points, x):
    s = cubic_spline_coeff(x_points, y_points)
    for i in range(len(x_points) - 1):
        if x_points[i] < x < x_points[i + 1]:
            result = (6 * s[i][0] * x) + (2 * s[i][1])
            return result

def cubic_spline_second_derivative(x_points, y_points):
    return lambda x: cubic_spline_second_derivative_result(x_points, y_points, x)



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


def first_order_spline_coeffs(x_points, y_points):
    h = [y_points[i + 1] - y_points[i] for i in range(len(y_points) - 1)]
    k = [x_points[i + 1] - x_points[i] for i in range(len(y_points) - 1)]
    a = [h[i] / k[i] for i in range(len(y_points) - 1)]
    b = [y_points[i] - x_points[i] * (h[i] / k[i]) for i in range(len(y_points) - 1)]
    return a, b

def first_order_spline_result(x_points, y_points, x):
    a, b = first_order_spline_coeffs(x_points, y_points)
    for i in range(len(x_points) - 1):
        if x_points[i] < x < x_points[i + 1]:
            return a[i] * x + b[i]

def first_order_spline(x_points, y_points):
    return lambda x: first_order_spline_result(x_points, y_points, x)

def hermite_result(x_points, y_points, tension, x):
    v_in = [tension * (y_points[i + 1] - y_points[i - 1]) for i in range(1, len(y_points) - 1)]
    v = [(y_points[1] - y_points[0]) / 2] + v_in + [(y_points[-1] - y_points[-2]) / 2]
    for i in range(len(x_points) - 1):
        if x_points[i] <= x <= x_points[i + 1]:
            new_x = (x - x_points[i]) / (x_points[i + 1] - x_points[i])
            H_0 = (2 * new_x ** 3) - (3 * new_x ** 2) + 1
            H_1 = (new_x ** 3) - (2 * new_x ** 2) + new_x
            H_2 = (-2 * new_x ** 3) + (3 * new_x ** 2)
            H_3 = (new_x ** 3) - (new_x ** 2)
            return H_0 * y_points[i] + H_1 * v[i] + H_2 * y_points[i + 1] + H_3 * v[i + 1]

def Hermite(x_points, y_points, tension):
    return lambda x: hermite_result(x_points, y_points, tension, x)

def step_function(x):
    if 0 <= x < 2.5:
        return 1
    if 2.5 <= x <=5:
        return -1
    return 0

x_points_1 = [0, 2, 3, 4, 7, 8]
y_points_1 = [4, 2, 8, 10, 4, -2]
x_points_2 = [3, 2, 2.5, 4, 5, 4]
y_points_2 = [4, 3, 1, 2, 3.5, 4.5]
x_points_3 = [0, 1, 2, 3, 4, 5]
y_points_3 = [1, 1, 1, -1, -1, -1]

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

