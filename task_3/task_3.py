import numpy as np
from tools import multiplying_matrix, LU_Decomposition, gauss, partial_pivoting, solve_with_LU, inverse_matrix, gauss_seidel

matrix = [[3, -3, 2, -4], [-2, -1, 3, -1], [5, -2, -3, 2], [-2, 4, 1, 2]]
vector = [7.9, -12.5, 18, -8.1]
A = [[4, 8, 4, 0], [1, 4, 7, 2], [1, 5, 4, -3], [1, 3, 0, -2]]


if __name__ == '__main__':
    print(f"Using Gauss Elimination the result is {gauss(matrix, vector)}")
    print(f"Using LU Decomposition the result is {solve_with_LU(matrix, vector)}")
    print(f"Using Gauss-Seidel the result is {gauss_seidel(matrix, vector)}")
    print(inverse_matrix(A))
    print(LU_Decomposition(A))
    A_1 = inverse_matrix(A)
    print(multiplying_matrix(A, A_1))