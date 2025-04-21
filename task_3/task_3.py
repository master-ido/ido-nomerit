def gauss(matrix, vector):
    # first line
    for j in range(1, len(matrix[0])):
        matrix[0][j] = matrix[0][j] / matrix[0][0]
    vector[0] = vector[0] / matrix[0][0]
    matrix[0][0] = 1
    # other lines
    for k in range(1, len(matrix)):
        for j in range(len(matrix)):
            matrix[k][j] = matrix[k][j] - (matrix[k][j] * matrix[k - 1][j])
        vector[k] = vector[k] - (matrix[k][j] * vector[k - 1])
    # finding the result
    result = [0 for i in range(len(vector))]
    result[-1] = vector[-1]
    for x in range(len(vector) - 1, 0, -1):
        result[x] = vector[x] - sum(matrix[i][i + 1] * result[i + 1] for i in range(x, len(vector) - 1))
    return result

matrix = [[3, -3, 2, -4], [-2, -1, 3, -1], [5, -2, -3, 2], [-2, 4, 1, 2]]
vector = [7.9, -12.5, 18, -8.1]
print(gauss(matrix, vector))




# # # a = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
# # # a[1][1] = a[1][1] / 5
# # # print(a)
#
a = "kop"
list = [i for i in range(len(a))]
# print(list)


