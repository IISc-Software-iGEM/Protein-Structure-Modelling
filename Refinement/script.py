import numpy as np

# Function to read the matrix from a file
def read_matrix_from_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        matrix = [list(map(float, line.split())) for line in lines]
    return matrix

# Function to transpose the matrix
def transpose_matrix(matrix):
    matrix_np = np.array(matrix)
    transpose_np = matrix_np.T
    return transpose_np.tolist()

# File path
file_path = 'test.txt'

# Read the matrix from the file
matrix = read_matrix_from_file(file_path)

# Transpose the matrix
transpose = transpose_matrix(matrix)

# Print the transposed matrix
for row in transpose:
    print(row[0], row[1], row[2])
