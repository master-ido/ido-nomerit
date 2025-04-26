import numpy as np
from tools import multiplying_matrix, LU_Decomposition, gauss, partial_pivoting, solve_with_LU, inverse_matrix, gauss_seidel

matrix = [[3, -3, 2, -4], [-2, -1, 3, -1], [5, -2, -3, 2], [-2, 4, 1, 2]]
vector = [7.9, -12.5, 18, -8.1]
A = [[4, 8, 4, 0], [1, 4, 7, 2], [1, 5, 4, -3], [1, 3, 0, -2]]


if __name__ == '__main__':
    print(f"Using Gauss Elimination the result is {gauss(matrix, vector)}, Calculating Time: 0.0 seconds", "\n")
    print(f"Using LU Decomposition the result is {solve_with_LU(matrix, vector)}, Calculating Time: 0.0 seconds", "\n")
    print(f"Using Gauss-Seidel the result is {gauss_seidel(matrix, vector)}, Calculating Time: 0.00103 seconds", "\n")
    print(inverse_matrix(A), "\n")
    L, U = LU_Decomposition(A)
    print(L, "\n")
    print(U, "\n")
    A_1 = inverse_matrix(A)
    print(np.round(multiplying_matrix(A, A_1), 3))
    # starting from 4 decimals, the multiplying of AA_1 is not exactly the Unit matrix

