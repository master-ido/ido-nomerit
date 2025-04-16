import math
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
    for i in np.arange(start, end, step):
        if bisection(function, i, i + step) == None:
            print(f"There are no roots between {i}, {i + step}")
        else:
            print(f"The root between {i}, {i + step} is {bisection(function, i, i + step)}")

def numeric_derivative(function, x):
    return (function(x + 10 ** -7) - function(x)) / (10 ** -7)


def newton_raphson(function, x):
    x_new = x - (function(x) / numeric_derivative(function, x))
    if abs(x_new - x) < 10 ** -7:
        return x_new
    return newton_raphson(function, x_new)

def find_roots_with_newton_raphson(function, start, end, step):
    for i in np.arange(start, end, step):
        print(f'Using the initial guess {i}, the root is {newton_raphson(function, i)}')


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
        if round(synthetic_division(b_n, i), 7) in roots:
            continue
        roots.append(round(synthetic_division(b_n, i), 7))
    return roots