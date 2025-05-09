#include <iostream>      // For input/output stream operations
#include <vector>        // To use std::vector for matrix representation
#include <cstdlib>       // For random number generation: rand(), srand()
#include <ctime>         // For seeding rand() with time()
#include <chrono>        // For measuring execution time
#include <unordered_set> // To avoid repeated column indices in sparse matrix
#include <tuple>         // Not used in current code, but included

using namespace std;

// ------------------------ DENSE ----------------------------------
// Function to perform dense matrix multiplication (standard triple-loop method)
void runDenseMultiplication(int n) {
    // Declare square matrices A, B and result matrix C initialized to 0
    vector<vector<int>> A(n, vector<int>(n));
    vector<vector<int>> B(n, vector<int>(n));
    vector<vector<int>> C(n, vector<int>(n, 0));

    // Seed the random generator and fill matrices A and B with values in [0, 9]
    srand(time(NULL));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            A[i][j] = rand() % 10;
            B[i][j] = rand() % 10;
        }

    // Start timing
    auto start = chrono::high_resolution_clock::now();

    // Perform matrix multiplication: C = A × B
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                C[i][j] += A[i][k] * B[k][j];

    // Stop timing
    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double>(end - start).count();

    // Print results
    cout << "Dense matrix multiplication completed.\n";
    cout << "Time taken: " << elapsed << " seconds.\n";
}

// ------------------------ SPARSE (CSR) ----------------------------

// Data structure representing a sparse matrix using CSR format
struct CSRMatrix {
    int rows, cols;               // Number of rows and columns in the matrix
    vector<int> values;           // Non-zero values in the matrix
    vector<int> col_indices;      // Column indices for each value
    vector<int> row_pointers;     // Row pointer array indicating where each row starts in 'values'
};

// Function to generate a random sparse matrix in CSR format
CSRMatrix generateSparseMatrix(int rows, int cols, double sparsity) {
    CSRMatrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.row_pointers.push_back(0); // The first row starts at index 0
    srand(time(NULL));             // Seed random number generator

    // Loop through each row
    for (int i = 0; i < rows; ++i) {
        unordered_set<int> used_cols; // To ensure no duplicate columns in the same row
        int non_zeros = static_cast<int>(cols * sparsity + 0.5); // Number of non-zero elements in this row

        // Generate random column indices and values
        for (int j = 0; j < non_zeros; ++j) {
            int col;
            do {
                col = rand() % cols; // Choose a random column
            } while (used_cols.count(col)); // Repeat until column is unique in this row
            used_cols.insert(col);

            mat.values.push_back(rand() % 10); // Random value between 0 and 9
            mat.col_indices.push_back(col);   // Store its column
        }

        // Store the current size of 'values' as the end of the row
        mat.row_pointers.push_back(mat.values.size());
    }

    return mat;
}

// Function to multiply two CSR-format sparse matrices: A × B = C
CSRMatrix multiplySparse(const CSRMatrix& A, const CSRMatrix& B) {
    // Ensure matrices are compatible for multiplication
    if (A.cols != B.rows) {  // check if A's columns match B's rows
        cerr << "Matrix size mismatch!\n";
        exit(1);
    }

    int C_rows = A.rows, C_cols = B.cols;

    // Create a temporary dense result matrix for computation for simplicity
    vector<vector<int>> temp(C_rows, vector<int>(C_cols, 0)); 

    // Multiply A and B using their sparse representation
    for (int i = 0; i < A.rows; ++i) {
        // Traverse through non-zero elements of row i in matrix A from row pointers
        for (int j = A.row_pointers[i]; j < A.row_pointers[i + 1]; ++j) {
            int A_col = A.col_indices[j]; // Column of the current A element
            int A_val = A.values[j];      // Value of the current A element

            // Traverse through non-zero elements of row A_col in matrix B
            for (int k = B.row_pointers[A_col]; k < B.row_pointers[A_col + 1]; ++k) {
                int B_col = B.col_indices[k]; // Column of the current B element
                int B_val = B.values[k];      // Value of the current B element

                // Multiply and accumulate in the temp matrix
                temp[i][B_col] += A_val * B_val;
            }
        }
    }

    // Now convert the temp matrix back to CSR format
    CSRMatrix C;
    C.rows = C_rows;
    C.cols = C_cols;
    C.row_pointers.push_back(0); // First row starts at index 0

    // Traverse the temp result and fill CSR structure
    for (int i = 0; i < C_rows; ++i) {
        for (int j = 0; j < C_cols; ++j) {
            if (temp[i][j] != 0) {
                C.values.push_back(temp[i][j]);
                C.col_indices.push_back(j);
            }
        }
        C.row_pointers.push_back(C.values.size()); // End index of current row
    }

    return C;
}

// Function to execute sparse multiplication and measure its performance
void runSparseMultiplication(int rows, int cols, double sparsity) {
    // Generate two compatible sparse matrices
    CSRMatrix A = generateSparseMatrix(rows, cols, sparsity);
    CSRMatrix B = generateSparseMatrix(cols, rows, sparsity); // Transposed dimensions for multiplication

    // Measure time for sparse matrix multiplication
    auto start = chrono::high_resolution_clock::now();
    CSRMatrix C = multiplySparse(A, B);
    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double>(end - start).count();

    cout << "Sparse matrix (CSR) multiplication completed.\n";
    cout << "Time taken: " << elapsed << " seconds.\n";
}

// ------------------------ MAIN ----------------------------------
// Main function to allow user to choose between dense and sparse multiplication
int main() {
    char choice;
    cout << "Matrix Multiplication Program\n";
    cout << "Choose matrix type:\n";
    cout << "[d] Dense\n[s] Sparse\n> ";
    cin >> choice;

    // User chooses dense multiplication
    if (choice == 'd' || choice == 'D') {
        int n;
        cout << "Enter size of square matrix (e.g., 3000): ";
        cin >> n;
        runDenseMultiplication(n);
    }
    // User chooses sparse multiplication
    else if (choice == 's' || choice == 'S') {
        int n;
        double sparsity;
        cout << "Enter size of square matrix (e.g., 3000): ";
        cin >> n;
        cout << "Enter sparsity (0.0 to 1.0): ";
        cin >> sparsity;
        runSparseMultiplication(n, n, sparsity);
    }
    // Invalid input
    else {
        cout << "Invalid choice. Use 'd' or 's'.\n";
    }

    return 0;
}