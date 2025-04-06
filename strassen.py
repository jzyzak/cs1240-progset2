import sys
#import time
#import numpy as np
#import matplotlib.pyplot as plt
import random
#from math import comb
import math

def next_power_of_two(n):
    if n == 0:
        return 1
    return 2 ** math.ceil(math.log2(n))


# Function for padding matrix when n isn't even
def pad_matrix(matrix, new_size):
    # Old padding method (did it for odd numbers but something was messed up)
    """padded = []

    # Adds zero column to end of matrix
    for row in matrix:
        padded.append(row + [0]) 

    # Adds zero row to end of matrix
    padded.append([0] * (n + 1)) 

    return padded"""
    # New padding method (just pads to the next power of 2 to minimize room for error when it comes to padding recursively like in old approach)
    old_n = len(matrix)
    padded = [[0]*new_size for _ in range(new_size)]
    for i in range(old_n):
        for j in range(old_n):
            padded[i][j] = matrix[i][j]
    return padded


# Removes padding so we can go back to original size matrix
def remove_padding(matrix, original_size):
    """
    n = len(matrix)
    for i in range(n):
        matrix[i] = matrix[i][:original_size]
    """

    return [row[:original_size] for row in matrix[:original_size]]


# Function for adding matrices
def add_matrix(A, B):
    n = len(A)
    result = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            result[i][j] = A[i][j] + B[i][j]
    return result


# Function for subtracting matrices (just the above function with "-" instead of "+")
def subtract_matrix(A,B):
    n = len(A)
    diff = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            diff[i][j] = A[i][j] - B[i][j]
    return diff


# Splits matrix into the 4 submatrices
def split_matrix(matrix):
    n = len(matrix)
    mid = n // 2

    A11 = [row[:mid] for row in matrix[:mid]]
    A12 = [row[mid:] for row in matrix[:mid]]
    A21 = [row[:mid] for row in matrix[mid:]]
    A22 = [row[mid:] for row in matrix[mid:]]
    return A11, A12, A21, A22


# Combines submatrices into one full matrix
def combining_results(AE_BG, AF_BH, CE_DG, CF_DH):
    n = len(AE_BG)
    #C = [[0]*n for _ in range(n)]
    C = [[0] * (2 * n) for _ in range(2 * n)]
    for i in range(n):
        for j in range(n):
            C[i][j] = AE_BG[i][j]
            C[i][j+n] = AF_BH[i][j]
            C[i+n][j] = CE_DG[i][j]
            C[i+n][j+n] = CF_DH[i][j]
    return C
    

# Standard matrix multiplication function
def standard_matrix_multiplication(A, B):
    n = len(A)
    C = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                C[i][j] += A[i][k] * B[k][j]
    return C


# Strassen's algorithm :)
def strassen(A, B, n0):

    n = len(A)
    next_n = next_power_of_two(n)
    
    # Pads to the next power of 2
    if next_n != n:
        A = pad_matrix(A, next_n)
        B = pad_matrix(B, next_n)

    # Handles the recursive part of Strassen's algorithm
    def strassen_recursion(A, B):
        n = len(A)
        if n <= n0 or n == 1:
            if n == 1:
                return [[A[0][0] * B[0][0]]]
            return standard_matrix_multiplication(A, B)

        A1, B1, C, D = split_matrix(A)
        E, F, G, H = split_matrix(B)

        P1 = strassen_recursion(A1, subtract_matrix(F, H))
        P2 = strassen_recursion(add_matrix(A1, B1), H)
        P3 = strassen_recursion(add_matrix(C, D), E)
        P4 = strassen_recursion(D, subtract_matrix(G, E))
        P5 = strassen_recursion(add_matrix(A1, D), add_matrix(E, H))
        P6 = strassen_recursion(subtract_matrix(B1, D), add_matrix(G, H))
        P7 = strassen_recursion(subtract_matrix(C, A1), add_matrix(E, F))

        AE_BG = add_matrix(subtract_matrix(add_matrix(P5, P4), P2), P6)
        AF_BH = add_matrix(P1, P2)
        CE_DG = add_matrix(P3, P4)
        CF_DH = add_matrix(subtract_matrix(add_matrix(P1, P5), P3), P7)

        return combining_results(AE_BG, AF_BH, CE_DG, CF_DH)

    C = strassen_recursion(A, B)

    if next_n != n:
        C = remove_padding(C, n)

    return C


# Function for counting triangles in a graph using Strassen's to get A^3
def triangle_count(A):
    A2 = strassen(A,A,16)
    A3= strassen(A2,A,16)

    total_triangles = sum(A3[i][i] for i in range(len(A3))) // 6
    return total_triangles


# Function for generating a random matrix of size n x n
def generate_matrix(n):
    return [[random.randint(0, 10) for _ in range(n)] for _ in range(n)]

"""
# Function for finding the crossover point between Strassen's and standard multiplication experimentally
def test_experimental_crossover(matrix_size=128):
    print(f"Testing experimental crossover for matrix size {matrix_size}×{matrix_size}")
    A = generate_matrix(matrix_size)
    B = generate_matrix(matrix_size)

    # Timing for standard matrix multiplication (just need to do once; explanation in write up)
    t0 = time.perf_counter()
    standard_matrix_multiplication(A, B)
    t1 = time.perf_counter()
    standard_time = t1 - t0
    print(f"Standard Time: {standard_time:.5f} s")

    # Testing Strassen's algorithm with different n0 values to try and find crossover point
    n0_vals = list(range(2, 50, 2))
    strassen_times = []

    for n0 in n0_vals:
        t0 = time.perf_counter()
        strassen(A, B, n0)
        t1 = time.perf_counter()
        duration = t1 - t0
        strassen_times.append(duration)
        print(f"n0 = {n0:2d} | Strassen Time: {duration:.5f} s")

    # Plot the results
    plt.figure(figsize=(10, 6))
    plt.plot(n0_vals, strassen_times, marker='o', label="Strassen")
    plt.axhline(standard_time, color='red', linestyle='--', label="Standard")
    plt.xlabel("Crossover Point n0")
    plt.ylabel("Runtime (seconds)")
    plt.title(f"Strassen vs Standard Multiplication (n = {matrix_size})")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Find experimental crossover
    for n0, t in zip(n0_vals, strassen_times):
        if t < standard_time:
            print(f"\n Experimental crossover point is around n0 = {n0}")
            return n0

    print("\n Strassen was never faster than standard for this size.")
    return None
"""

# Function for generating the random graph based on edge probabilities
def generate_random_graph(n, p):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() <= p:
                A[i][j] = 1
                A[j][i] = 1
    return A

"""
# Function for counting triangles on the random graphs with 1024 vertices with varying edge probabilities
def run_triangle_experiment(n=1024, ps=[0.01, 0.02, 0.03, 0.04, 0.05]):
    actual_counts = []
    expected_counts = []

    for p in ps:
        print(f"Running for p = {p:.2f}")
        A = generate_random_graph(n, p)
        actual = triangle_count(A)
        expected = comb(n, 3) * (p ** 3)
        actual_counts.append(actual)
        expected_counts.append(expected)
        print(f"Actual: {actual}, Expected: {expected:.2f}")

    plt.figure(figsize=(10, 6))
    plt.plot(ps, actual_counts, marker='o', label="Actual Triangle Count")
    plt.plot(ps, expected_counts, marker='x', linestyle='--', label="Expected Triangle Count")
    plt.xlabel("Edge Probability (p)")
    plt.ylabel("Triangle Count")
    plt.title("Actual vs Expected Triangle Counts in Random Graphs (n = 1024)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
"""

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 strassen.py <flag> <dimension> <inputfile>")
        sys.exit(1)

    flag = int(sys.argv[1])
    d = int(sys.argv[2])
    inputfile = sys.argv[3]

    with open(inputfile, 'r') as f:
        nums = [int(line.strip()) for line in f.readlines()]

    A_flat = nums[:d*d]
    B_flat = nums[d*d:]
    A = [A_flat[i*d:(i+1)*d] for i in range(d)]
    B = [B_flat[i*d:(i+1)*d] for i in range(d)]

    C = strassen(A, B, n0=16)
    for i in range(d):
        print(C[i][i])
    #test_experimental_crossover(matrix_size=128)
    #A = [[1,1,1],[1,0,1],[0,0,1]]
    #B = [[1,0,1],[0,0,1],[1,1,1]]
    #print(strassen(A,B,0))
    #run_triangle_experiment(n=1024, ps=[0.01, 0.02, 0.03, 0.04, 0.05])