import sys
import time
import numpy as np

# Testing values for the theoretical cost

def strassen_cost(n,n0,c=1,c_prime=1):
    leaf_cost = (n/n0)**np.log2(7) * c_prime*(n0**3)
    internal_cost = (4/3)*c*(n**2) * ((n/n0)**np.log2(7/4) - 1)
    return leaf_cost + internal_cost

# Pads matrix to give it even dimensions/handle odd dimension case
def pad_matrix(matrix,n):
    for row in matrix:
        row.append(0)
    
    matrix.append([0]*(n+1))
    return matrix

def remove_padding(matrix, original_size):
    n = len(matrix)
    for i in range(n):
        matrix[i] = matrix[i][:original_size]

    return matrix


# Splits matrix into the 4 submatrices
def split_matrix(matrix):
    n = len(matrix)
    mid = n // 2

    A11 = [row[:mid] for row in matrix[:mid]]
    A12 = [row[mid:] for row in matrix[:mid]]
    A21 = [row[:mid] for row in matrix[mid:]]
    A22 = [row[mid:] for row in matrix[mid:]]
    return A11, A12, A21, A22

# Function for adding matrices
def add_matrix(A,B):
    n = len(A)
    sum = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            sum[i][j] = A[i][j] + B[i][j]
    return sum

def combining_results(AE_BG, AF_BH, CE_DG, CF_DH):
    
    pass

# Function for subtracting matrices (just the above function with "-" instead of "+")
def subtract_matrix(A,B):
    n = len(A)
    diff = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            sum[i][j] = A[i][j] - B[i][j]
    return diff


# Standard matrix multiplication function
def standard_matrix_multiplication(A, B):
    n = A.shape[0]

    C = np.zeros((n, n))
    for i in range (0, n):
        for j in range (0, n):
            C[i][j] = 0
            for k in range (0, n):
                C[i][j] += A[i][k] * B[k][j]
    pass

def strassen(A,B, n0):
    n = A.shape[0]
    padded = False

    # Base Case for when standard matrix multiplication is better
    if n <= n0:
        return standard_matrix_multiplication(A, B)
    
    # Adds empty padding for when n is not even
    if n0%2==1:
        padded = True
        A = pad_matrix(A, n)
        B = pad_matrix(B, n)
        n += 1

    # Split matrices intro components
    A_split, B_split, C, D = split_matrix(A)
    E, F, G, H = split_matrix(B)

    P1 = strassen(A_split, subtract_matrix(F, H), n0)
    P2 = strassen(add_matrix(A_split, B_split), H, n0)
    P3 = strassen(add_matrix(C, D), E, n0)
    P4 = strassen(D, subtract_matrix(G, E), n0)
    P5 = strassen(add_matrix(A_split, D), add_matrix(E, H), n0)
    P6 = strassen(subtract_matrix(B_split, D), add_matrix(G, H), n0)
    P7 = strassen(subtract_matrix(C, A_split), add_matrix(E, F), n0)
    
    AE_BG = add_matrix(add_matrix(-P2, P4),add_matrix(P5,P6))
    AF_BH = add_matrix(P1, P2)
    CE_DG = add_matrix(P3, P4)
    CF_DH = add_matrix(subtract_matrix(P1, P3), add_matrix(P5,P7))
    pass

n=1024
n0_values=[2**i for i in range(1,10)]

for n0 in n0_values:
    total_cost = strassen_cost(n,n0)
    print(f"n0 = {n0:3d} -> total_cost = {total_cost}")