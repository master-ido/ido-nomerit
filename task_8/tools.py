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


def present_function(function, start, end):
    x = np.linspace(start, end, 500)
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
    v = [(y_points[1] - y_points[0]) * tension] + v_in + [(y_points[-1] - y_points[-2]) * tension]
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

def trapeze_method(function, limits, dividing_steps):
    x = [limits[0] + ((limits[1] - limits[0]) / (dividing_steps - 1)) * i for i in range(dividing_steps)]
    h = (limits[1] - limits[0]) / (dividing_steps - 1)
    return sum((h / 2) * (function(x[i]) + function(x[i + 1])) for i in range(len(x) - 1))

def simpson_1_3_combined(function, limits, dividing_steps):
    x = np.linspace(limits[0], limits[1], dividing_steps)
    h = (limits[1] - limits[0]) / dividing_steps
    result = function(x[0]) + function(x[-1])
    for i in range(1, len(x) - 1):
        if i % 2 == 0:
            result += (function(x[i]) * 2)
        else:
            result += (function(x[i]) * 4)
    return (h / 3) * result

def simpson_3_8_combined(function, limits, analytic_result, iteration=1):
    x = np.linspace(limits[0], limits[1], iteration * 3 + 1)
    h = (limits[1] - limits[0]) / (iteration * 3)
    result = function(x[0]) + function(x[-1])
    for i in range(1, len(x) - 1):
        if i % 3 == 0:
            result += (function(x[i]) * 2)
        else:
            result += (function(x[i]) * 3)
    if abs(((3 * h) / 8) * result - analytic_result) > 10 ** -5:
        return simpson_3_8_combined(function, limits, analytic_result, iteration + 1)
    return ((3 * h) / 8) * result, iteration

def trapeze_combined(function, limits, n):
    h = (limits[1] - limits[0]) / n
    x = np.linspace(limits[0], limits[1], n)
    result = function(limits[0]) + function(limits[1]) + sum(2 * function(x[i]) for i in range(1, len(x) - 1))
    return (h / 2) * result

def romberg_loop(function, limits, analytic_result):
    mat = [[]]
    for i in range(200):
        mat.append([0] * (len(mat[0]) + 1))
        for row in mat[:-1]:
            row.append(0)
        n = 2 ** i
        mat[i][0] = trapeze_combined(function, limits, n)
        for k in range(1, i + 1):
            mat[i][k] = (4 ** k * mat[i][k - 1] - mat[i - 1][k - 1]) / ((4 ** k) - 1)
        if i > 0 and abs(mat[i][i] - analytic_result) < 10 ** -3:
            return mat[i][i], i

def romberg_recursion(function, limits, analytic_result, mat=None, i=0):
    if mat is None:
        mat = [[]]
    mat.append([0] * (len(mat[0]) + 1))
    for row in mat[:-1]:
        row.append(0)
    n = 2 ** i
    mat[i][0] = trapeze_combined(function, limits, n)
    for k in range(1, i + 1):
        mat[i][k] = (4 ** k * mat[i][k - 1] - mat[i - 1][k - 1]) / ((4 ** k) - 1)
    if abs(mat[i][i] - analytic_result) > 10 ** -3:
        return romberg_recursion(function, limits, analytic_result, mat, i + 1)
    return mat[i][i], i

def quad_function(function, limits):
    a = ((limits[1] - limits[0]) / 2)
    b = ((limits[1] + limits[0]) / 2)
    return lambda t: (function((a * t) + b)) * a

quad_table = [[0.57735026, -0.57735026, 1, 1], [0.77459667, -0.77459667, 0, 0.88888889, 0.55555555, 0.55555555], [0.86113631, -0.86113631, 0.33998104, -0.33998104, 0.65214515, 0.65214515, 0.34785485, 0.34785485], [0.90617985, -0.90617985, 0.53846931, -0.53846931, 0, 0.56888889, 0.47862867, 0.47862867, 0.23692689, 0.23692689], [0.93246951, -0.93246951, 0.66120939, -0.66120939, 0.23861918, -0.23861918, 0.46791393, 0.46791393, 0.36076157, 0.36076157, 0.17132449, 0.17132449], [0.94910791, -0.94910791, 0.74153119, -0.74153119, 0.40584515, -0.40584515, 0, 0.41795918, 0.38183005, 0.38183005, 0.27970539, 0.27970539, 0.12948497, 0.12948497], [0.96028986, -0.96028986, 0.79666648, -0.79666648, 0.52553241, -0.52553241, 0.18343464, -0.18343464, 0.36268378, 0.36268378, 0.31370665, 0.31370665, 0.22238103, 0.22238103, 0.10122854, 0.10122854], [0.97390653, -0.97390653, 0.86506337, -0.86506337, 0.67940957, -0.67940957, 0.43339539, -0.43339539, 0.14887434, -0.14887434, 0.29552422, 0.29552422, 0.26926672, 0.26926672, 0.21908636, 0.21908636, 0.14945135, 0.14945135, 0.06667134, 0.06667134]]

def open_quad(function, limits, analytic_result, iteration = 0):
    table = quad_table
    g = quad_function(function, limits)
    result = sum(g(table[iteration][i]) * table[iteration][len(table[iteration]) - 1 - i] for i in range(len(table[iteration]) // 2))
    if abs(result - analytic_result) > 10 ** -5 and iteration < 7:
        return open_quad(function, limits, analytic_result, iteration + 1)
    return result, iteration

def simpson_3_over_8(function, limits):
    h = (limits[1] - limits[0]) / 6
    a_1 = (limits[0] + limits[1]) / 2
    return h * (function(limits[0]) + function(a_1) + function(limits[1]))

def finite_diff_order_4_error(function, x, h):
    forward = (-25 * function(x) + 48 * function(x + h) - 36 * function(x + (2 * h)) + 16 * function(x + (3 * h)) - 3 * function(x + (4 * h))) / (12 * h)
    backward = (25 * function(x) - 48 * function(x - h) + 36 * function(x - (2 * h)) - 16 * function(x - (3 * h)) + 3 * function(x - (4 * h))) / (12 * h)
    centered = (8 * function(x + h) - 8 * function(x - h) - function(x + (2 * h)) + function(x - (2 * h))) / (12 * h)
    return forward, backward, centered

# def euler_method(dy_dx, initial_x, initial_y, final_x, h, results=None):
#     if results is None:
#         results = [(initial_x, initial_y)]
#     results.append((initial_x + h, initial_y + h * dy_dx(initial_x, initial_y)))
#     if round(initial_x + h, 3) == final_x:
#         return results
#     return euler_method(dy_dx, initial_x + h, initial_y + h * dy_dx(initial_x, initial_y), final_x, h, results)
def euler_method(dy_dx, initial_x, initial_y, final_x, h, x_results=None, y_results=None):
    if x_results is None:
        x_results = [initial_x]
        y_results = [initial_y]
    x_results.append(initial_x + h)
    y_results.append(initial_y + h * dy_dx(initial_x, initial_y))
    if round(initial_x + h, 3) == final_x:
        return x_results, y_results
    return euler_method(dy_dx, initial_x + h, initial_y + h * dy_dx(initial_x, initial_y), final_x, h, x_results, y_results)

def euler_2nd_order(dv_dt, initial_t, initial_v, initial_y, final_t, h, t_results=None, y_results=None, v_results=None):
    if t_results is None:
        t_results = [initial_t]
        v_results = [initial_v]
        y_results = [initial_y]
    t_results.append(initial_t + h)
    v_results.append(initial_v + h * dv_dt(initial_t, initial_v, initial_y))
    y_results.append(initial_y + h * initial_v)
    if round(initial_t + h, 3) == final_t:
        return [round(t, 4) for t in t_results], [round(v, 4) for v in v_results], [round(y, 4) for y in y_results]
    return euler_2nd_order(dv_dt, initial_t + h, initial_v + h * dv_dt(initial_t, initial_v, initial_y), initial_y + h * dy_dt(initial_v), final_t, h, t_results, y_results, v_results)

def present_euler_2nd(dv_dt, initial_t, initial_v, initial_y, final_t, h):
    t_vals, v_vals, y_vals = euler_2nd_order(dv_dt, initial_t, initial_v, initial_y, final_t, h)
    plt.plot(t_vals, y_vals, marker='o', linestyle='-', color='blue', label="Euler's 2nd Method")
    plt.xlabel('t')
    plt.ylabel('y')
    plt.title("Euler's 2nd Method Approximation")
    plt.grid(True)
    plt.legend()
    plt.show()
def dv_dt(t, v, y):
    return -20 * y
def dy_dt(v):
    return v
# print(euler_2nd_order(dv_dt, 0, 0, 0.7, 2.5, 0.1))
# present_euler_2nd(dv_dt, 0, 0, 0.7, 2.5, 0.01)
def dm_dx(m,x):
    return 10 - 2 * x
# print(euler_method(dm_dx, 0, 0, 1, 0.05))

def runga_kutta_2nd_order(dy_dx, x, y, final_x, h, x_results=None, y_results=None):
    if x_results is None:
        x_results = [x]
        y_results = [y]
    f_i = dy_dx(x, y)
    y_i = (h / 2) * (f_i + dy_dx(x + h, y + h * f_i)) + y
    x_results.append(round(x + h, 5))
    y_results.append(y_i)
    if round(x + h, 3) == final_x:
        return x_results, y_results
    return runga_kutta_2nd_order(dy_dx, x + h, y_i, final_x, h, x_results, y_results)

# print(runga_kutta_2_order(dm_dx, 0, 0, 1, 0.05))

def runga_kutta_4th_order(dy_dx, x, y, final_x, h, x_results=None, y_results=None):
    if x_results is None:
        x_results = [x]
        y_results = [y]
    F1 = dy_dx(x, y)
    F2 = dy_dx(x + (h / 2), y + (h / 2) * F1)
    F3 = dy_dx(x + (h / 2), y + (h / 2) * F2)
    F4 = dy_dx(x + h, y + h * F3)
    y_i = y + (h / 6) * (F1 + 2 * (F2 + F3) + F4)
    x_results.append(round(x + h, 5))
    y_results.append(y_i)
    if round(x + h, 3) == final_x:
        return x_results, y_results
    return runga_kutta_4th_order(dy_dx, x + h, y_i, final_x, h, x_results, y_results)


def runga_kutta_4_order_2nd(dv_dt, initial_t, initial_v, initial_y, final_t, h, t_results=None, v_results=None, y_results=None):
    if t_results is None:
        t_results = [initial_t]
        v_results = [initial_v]
        y_results = [initial_y]
    F1_y = initial_v
    F1_v = dv_dt(initial_t, initial_v, initial_y)
    F2_y = initial_v + (h / 2) * F1_v
    F2_v = dv_dt(initial_t + (h / 2), initial_v + (h / 2) * F1_v, initial_y + (h / 2) * F1_y)
    F3_y = initial_v + (h / 2) * F2_v
    F3_v = dv_dt(initial_t + (h / 2), initial_v + (h / 2) * F2_v, initial_y + (h / 2) * F2_y)
    F4_y = initial_v + h * F3_v
    F4_v = dv_dt(initial_t + h, initial_v + h * F3_v, initial_y + h * F3_y)
    y_next = initial_y + (h / 6) * (F1_y + 2 * (F2_y + F3_y) + F4_y)
    v_next = initial_v + (h / 6) * (F1_v + 2 * (F2_v + F3_v) + F4_v)
    t_results.append(round(initial_t + h, 5))
    v_results.append(v_next)
    y_results.append(y_next)
    if round(initial_t + h, 3) == final_t:
        return t_results, v_results, y_results
    return runga_kutta_4_order_2nd(dv_dt, initial_t + h, v_next, y_next, final_t, h, t_results, v_results, y_results)
# print(runga_kutta_4_order_2nd(dv_dt, 0, 0, 0.7, 2.5, 0.1))
def present_RK4_2nd(dv_dt, initial_t, initial_v, initial_y, final_t, h):
    t_vals, v_vals, y_vals = runga_kutta_4_order_2nd(dv_dt, initial_t, initial_v, initial_y, final_t, h)
    plt.plot(t_vals, y_vals, marker='o', linestyle='-', color='blue', label="RK4 2nd - order Method")
    plt.xlabel('t')
    plt.ylabel('y')
    plt.title("RK4 Method Approximation")
    plt.grid(True)
    plt.legend()
    plt.show()
# present_RK4_2nd(dv_dt, 0, 0, 0.7, 2.5, 0.1)

def present_euler(dy_dx, initial_x, initial_y, final_x, h):
    x_vals, y_vals = euler_method(dy_dx, initial_x, initial_y, final_x, h)
    plt.plot(x_vals, y_vals, marker='o', linestyle='-', color='blue', label="Euler's Method")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title("Euler's Method Approximation")
    plt.grid(True)
    plt.legend()
    plt.show()

def present_rk_2(dy_dx, initial_x, initial_y, final_x, h):
    x_vals, y_vals = runga_kutta_2nd_order(dy_dx, initial_x, initial_y, final_x, h)
    plt.plot(x_vals, y_vals, marker='o', linestyle='-', color='blue', label="RK2 Method")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title("RK2 Approximation")
    plt.grid(True)
    plt.legend()
    plt.show()

def present_rk_4(dy_dx, initial_x, initial_y, final_x, h):
    x_vals, y_vals = runga_kutta_4_order(dy_dx, initial_x, initial_y, final_x, h)
    plt.plot(x_vals, y_vals, marker='o', linestyle='-', color='blue', label="Euler's Method")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title("Euler's Method Approximation")
    plt.grid(True)
    plt.legend()
    plt.show()
# present_rk_4(dm_dx, 0,0,10,0.25)

# present_rk_2(dm_dx, 0,0,10,0.25)
# present_euler(dm_dx, 0, 0, 10, 0.25)

def adams_4_with_rk_4(dy_dx, x, y, final_x, h, x_results=None, y_results=None):
    if x_results is None:
        x_results, y_results = runga_kutta_4th_order(dy_dx, x, y, x + (3 * h), h)
    f1 = dy_dx(x_results[-1], y_results[-1])
    f2 = dy_dx(x_results[-2], y_results[-2])
    f3 = dy_dx(x_results[-3], y_results[-3])
    f4 = dy_dx(x_results[-4], y_results[-4])
    y_pre = y_results[-1] + (h / 24) * (55 * f1 - 59 * f2 + 37 * f3 - 9 * f4)
    y_cor_1 = y_results[-1] + (h / 24) * (9 * dy_dx(x_results[-1] + h, y_pre) + 19 * f1 - 5 * f2 + f3)
    y_cor_2 = y_results[-1] + (h / 24) * (9 * dy_dx(x_results[-1] + h, y_cor_1) + 19 * f1 - 5 * f2 + f3)
    y_cor_3 = y_results[-1] + (h / 24) * (9 * dy_dx(x_results[-1] + h, y_cor_2) + 19 * f1 - 5 * f2 + f3)
    x_results.append(x_results[-1] + h)
    y_results.append(y_cor_3)
    if round(x_results[-1], 3) == final_x:
        return x_results, y_results
    return adams_4_with_rk_4(dy_dx, x + h, y_cor_3, final_x, h, x_results, y_results)

def dy_dx(x, y):
    return (-2 * y) / (1 + x)

x_results, y_results = adams_4_with_rk_4(dy_dx, 0, 2, 10, 0.5)
# print(x_results, y_results)
# print(adams_4_with_rk_4(dy_dx, 0, 2, 10, 0.5))
def present_adams_4_with_rk_4(dy_dx, initial_x, initial_y, final_x, h):
    x_vals, y_vals = adams_4_with_rk_4(dy_dx, initial_x, initial_y, final_x, h)
    plt.plot(x_vals, y_vals, marker='o', linestyle='-', color='blue', label="Euler's Method")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title("Adam's Approximation")
    plt.grid(True)
    plt.legend()
    plt.show()
# present_adams_4_with_rk_4(dy_dx, 0, 2, 10, 0.5)
def leapfrog(dv_dt, initial_t, initial_v, initial_y, h, final_t, t_results=None, y_results=None):
    if t_results is None:
        t_results = [initial_t]
        y_results = [initial_y]
    half_step_v = initial_v + (h / 2) * dv_dt(initial_t, initial_v, initial_y)
    next_y = initial_y + half_step_v * h
    t_results.append(round(initial_t + h, 4))
    y_results.append(round(next_y, 5))
    if round(initial_t + h, 3) == final_t:
        return t_results, y_results
    full_step_v = half_step_v + (h / 2) * dv_dt(initial_t + h, half_step_v, next_y)
    return leapfrog(dv_dt, initial_t + h, full_step_v, next_y, h, final_t, t_results, y_results)
def d2x_dt2(y):
    return - (40 * y) / 2
initial_y = 0.7
initial_v = 0
# print(leapfrog(d2x_dt2, 0, 0, 0.7, 0.1, 2.5))

def present_leapfrog(dv_dt, initial_t, initial_v, initial_y, h, final_t):
    t_vals, y_vals = leapfrog(dv_dt, initial_t, initial_v, initial_y, h, final_t)
    plt.plot(t_vals, y_vals, marker='o', linestyle='-', color='blue', label="Leapfrog Method")
    plt.xlabel('t')
    plt.ylabel('y')
    plt.title("Leapfrog Method Approximation")
    plt.grid(True)
    plt.legend()
    plt.show()
# present_leapfrog(d2x_dt2, 0, 0, 0.7, 0.05, 2.5)

#
# def dy_dx(x, y):
#     return 4 * (x ** 2)
# initial_x = 1
# start_con = 1
# final_x = 1.1
# h = 0.05
# # print(euler_method(dy_dx, initial_x, start_con, final_x, h))
# print(runga_kutta_4_order(dy_dx, initial_x, start_con, final_x, h))
