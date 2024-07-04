import numpy as np
from scipy.sparse.linalg import eigs

def read_matrix(filename):
    """
    Reads a 24x24 matrix from a text file.
    Each row in the file should correspond to a row in the matrix.
    """
    return np.loadtxt(filename)

def compute_eigenvalues_and_vectors(G, dim):
    """
    Computes the largest 'dim' eigenvalues and eigenvectors of the matrix G.
    """
    eigenvalues, eigenvectors = eigs(G, k=dim, which='LM')
    
    # Ensure results are real (due to numerical errors, eigs might return small imaginary parts)
    eigenvalues = np.real(eigenvalues)
    eigenvectors = np.real(eigenvectors)
    
    return eigenvalues, eigenvectors

def main():
    filename = 'matrix.txt'  # Change to your filename
    G = read_matrix(filename)
    
    # Check if the matrix is 24x24
    if G.shape != (24, 24):
        print("Error: The matrix is not 24x24.")
        return
    
    dim = 3

    # Compute the eigenvalues and eigenvectors
    eigenvalues, eigenvectors = compute_eigenvalues_and_vectors(G, dim)

    print("Top 3 Eigenvalues:\n", eigenvalues)
    print("Top 3 Eigenvectors:\n", eigenvectors)

if __name__ == "__main__":
    main()
