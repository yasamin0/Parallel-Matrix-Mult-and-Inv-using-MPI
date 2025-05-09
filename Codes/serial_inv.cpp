#include <iostream>     // For input/output operations
#include <vector>       // To use 2D vector for matrix representation
#include <cstdlib>      // For rand(), srand()
#include <ctime>        // For time() used to seed srand()
#include <chrono>       // For measuring execution time

using namespace std;

// Function to compute the inverse of a square matrix using Gauss-Jordan elimination
void gaussJordanInversion(vector<vector<double>>& matrix) {
    int N = matrix.size();  // Matrix size (N x N)

    // Create an identity matrix of the same size
    vector<vector<double>> identity(N, vector<double>(N, 0.0));

    // Fill identity matrix (diagonal elements set to 1)
    for (int i = 0; i < N; i++) {
        identity[i][i] = 1.0;
    }

    // Begin Gauss-Jordan elimination process
    for (int i = 0; i < N; i++) {
        double pivot = matrix[i][i];  // Current pivot element , pivot onsore qotre alie dar satre i

        // Check if pivot is zero (singular matrix)
        if (pivot == 0.0) {
            cerr << "Matrix is singular or nearly singular.\n";
            return; // Cannot compute inverse
        }

        // Normalize the pivot row (divide entire row by pivot)
        for (int j = 0; j < N; j++) {
            matrix[i][j] /= pivot;
            identity[i][j] /= pivot;
        }

        // Eliminate the current column in all other rows
        for (int k = 0; k < N; k++) {
            if (k != i) { // Skip the pivot row
                double factor = matrix[k][i]; // Element to eliminate
                for (int j = 0; j < N; j++) {
                    matrix[k][j] -= factor * matrix[i][j];       // Make it zero in original matrix
                    identity[k][j] -= factor * identity[i][j];   // Apply same operation on identity matrix
                }
            }
        }
    }

    // After full elimination, identity becomes the inverse
    matrix = identity;
}
int main() {
    cout << "Matrix Inversion Program (Gauss-Jordan)\n";

    int n;
    cout << "Enter size of square matrix (e.g., 1000): ";
    cin >> n;

    // Declare n x n matrix and fill it with random values in [1, 10]
    vector<vector<double>> matrix(n, vector<double>(n));
    srand(time(NULL));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            matrix[i][j] = static_cast<double>(rand() % 10 + 1); // Avoid zero pivots

    // Record time before inversion
    auto start = chrono::high_resolution_clock::now();

    // Perform matrix inversion
    gaussJordanInversion(matrix);

    // Record time after inversion
    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double>(end - start).count();

    // Report completion and time taken
    cout << "Matrix inversion completed.\n";
    cout << "Time taken: " << elapsed << " seconds.\n";

    return 0;
}
