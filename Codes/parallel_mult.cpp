#include <iostream>         // For input/output
#include <vector>           // For using std::vector
#include <cstdlib>          // For rand(), srand()
#include <ctime>            // For seeding random generator
#include <unordered_set>    // For generating sparse matrix with unique columns
#include <mpi.h>            // For MPI (Message Passing Interface)

using namespace std;

// ------------------- Dense Generator -------------------
// Fills a matrix (stored in 1D vector) with random integers [0, 9]
void generateDenseMatrix(vector<int>& mat, int n) {
    srand(time(NULL));
    for (int i = 0; i < n * n; ++i)
        mat[i] = rand() % 10;
}

// ------------------- Sparse Generator -------------------
// Creates a sparse matrix by filling a limited number of non-zero values per row
// The matrix is stored as a 1D vector in row-major format
void generateSparseMatrixAsDense(vector<int>& mat, int n, double sparsity) {
    srand(time(NULL));
    fill(mat.begin(), mat.end(), 0); // Initialize entire matrix with 0

    int non_zeros_per_row = static_cast<int>(n * sparsity); // Number of non-zero elements per row
    for (int i = 0; i < n; ++i) {
        unordered_set<int> cols;
        while (cols.size() < non_zeros_per_row) {
            int col = rand() % n;
            mat[i * n + col] = rand() % 10; // Assign random value in allowed range
            cols.insert(col);               // Ensure column uniqueness
        }
    }
}

// ------------------- Multiply Block -------------------
// Performs matrix multiplication for a subset of rows (startRow to endRow)
// A and B are full matrices stored in 1D (row-major)
// Returns a submatrix of the result corresponding to the assigned rows
vector<int> multiplyBlock(const vector<int>& A, const vector<int>& B, int n, int startRow, int endRow) {
    vector<int> local_result((endRow - startRow) * n, 0); // Result for this block of rows

    for (int i = startRow; i < endRow; ++i) {
        for (int j = 0; j < n; ++j) {
            int sum = 0;
            for (int k = 0; k < n; ++k) {
                sum += A[i * n + k] * B[k * n + j]; // Standard matrix multiplication
            }
            local_result[(i - startRow) * n + j] = sum; // Store result relative to block
        }
    }
    return local_result;
}

// ------------------- Main Program -------------------
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv); // Initialize MPI environment

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get current process rank . shomare pardazande
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get total number of processes

    int n;             // Size of the square matrix
    char type;         // 'd' for dense, 's' for sparse
    double sparsity = 0.0;

    // Only the root process (rank 0) handles user input
    if (rank == 0) {
        cout << "Parallel Matrix Multiplication Program (MPI)\n";
        cout << "Choose matrix type:\n[d] Dense\n[s] Sparse\n> ";
        cin >> type;

        cout << "Enter size of square matrix (e.g., 3000): ";
        cin >> n;

        if (type == 's' || type == 'S') {
            cout << "Enter sparsity (0.0 to 1.0): ";
            cin >> sparsity;
        }
    }

    // Broadcast matrix type, size, and sparsity to all processes
    MPI_Bcast(&type, 1, MPI_CHAR, 0, MPI_COMM_WORLD); // (type , faqat 1 meqdar , noe dadeh , pardazande asli ke ettelaat pakhsh mikone , goruhe hame pardazandeha)
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sparsity, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Prepare matrices A and B
    vector<int> A(n * n, 0), B(n * n, 0);

    // Root process initializes the matrices . tolide matrices dar pardazande asli
    if (rank == 0) {
        if (type == 'd' || type == 'D') {
            generateDenseMatrix(A, n);
            generateDenseMatrix(B, n);
        } else {
            generateSparseMatrixAsDense(A, n, sparsity);
            generateSparseMatrixAsDense(B, n, sparsity);
        }
    }

    // Broadcast full matrices to all processes
    MPI_Bcast(A.data(), n * n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(B.data(), n * n, MPI_INT, 0, MPI_COMM_WORLD);

    // Divide the workload (rows) evenly among processes
    int rows_per_proc = n / size;
    int start_row = rank * rows_per_proc;
    int end_row = (rank == size - 1) ? n : start_row + rows_per_proc; // Handle remainder in last process

    MPI_Barrier(MPI_COMM_WORLD); // Sync before timing , hame pardazandeha montazer mimunan ta hame beresan be in noqte
    double start_time = MPI_Wtime(); // Start total timing

    // Each process performs local matrix multiplication
    vector<int> local_result = multiplyBlock(A, B, n, start_row, end_row);
    double computation_end = MPI_Wtime(); // Mark end of computation

    // Root process will collect all sub-results
    vector<int> full_result;
    if (rank == 0)
        full_result.resize(n * n);

    // Gather all sub-results from processes into the final result matrix
    MPI_Gather(local_result.data(), (end_row - start_row) * n, MPI_INT, //MPI_Gather(sendbuf, sendcount, sendtype,recvbuf, recvcount, recvtype,root, comm);
               full_result.data(), (end_row - start_row) * n, MPI_INT,
               0, MPI_COMM_WORLD); //shomare pardazande asli ke natijehaye hame pardazandeha ro jam mikone , kolle pardazandeha

    MPI_Barrier(MPI_COMM_WORLD); // Sync before finishing
    double end_time = MPI_Wtime(); // End total timing

    // Output results on root
    if (rank == 0) {
        double total_time = end_time - start_time;
        double comp_time = computation_end - start_time;
        double comm_time = total_time - comp_time;
        double comm_percent = (comm_time / total_time) * 100.0;

        if (type == 'd' || type == 'D')
            cout << "Dense matrix multiplication completed.\n";
        else
            cout << "Sparse matrix multiplication (via dense) completed.\n";

        cout << "Total Time: " << total_time << " seconds\n";
        cout << "Computation Time: " << comp_time << " seconds\n";
        cout << "Communication Time: " << comm_time << " seconds\n";
        cout << "Communication %: " << comm_percent << "%\n";
    }

    MPI_Finalize(); // Finalize MPI environment
    return 0;
}
