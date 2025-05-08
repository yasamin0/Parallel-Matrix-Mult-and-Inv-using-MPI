#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <chrono>

using namespace std;

void gaussJordanInversion(vector<vector<double>>& matrix) {
    int N = matrix.size();
    vector<vector<double>> identity(N, vector<double>(N, 0.0));

    // Initialize identity matrix
    for (int i = 0; i < N; i++) {
        identity[i][i] = 1.0;
    }

    // Perform Gauss-Jordan elimination
    for (int i = 0; i < N; i++) {
        double pivot = matrix[i][i];
        if (pivot == 0.0) {
            cerr << "Matrix is singular or nearly singular.\n";
            return;
        }

        // Normalize pivot row
        for (int j = 0; j < N; j++) {
            matrix[i][j] /= pivot;
            identity[i][j] /= pivot;
        }

        // Eliminate other rows
        for (int k = 0; k < N; k++) {
            if (k != i) {
                double factor = matrix[k][i];
                for (int j = 0; j < N; j++) {
                    matrix[k][j] -= factor * matrix[i][j];
                    identity[k][j] -= factor * identity[i][j];
                }
            }
        }
    }

    // Store the inverse
    matrix = identity;
}

int main() {
    cout << "🧮 Matrix Inversion Program (Gauss-Jordan)\n";
    int n;
    cout << "Enter size of square matrix (e.g., 1000): ";
    cin >> n;

    vector<vector<double>> matrix(n, vector<double>(n));
    srand(time(NULL));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            matrix[i][j] = static_cast<double>(rand() % 10 + 1); // Avoid zero pivots

    auto start = chrono::high_resolution_clock::now();
    gaussJordanInversion(matrix);
    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double>(end - start).count();

    cout << "✅ Matrix inversion completed.\n";
    cout << "⏱ Time taken: " << elapsed << " seconds.\n";

    return 0;
}
