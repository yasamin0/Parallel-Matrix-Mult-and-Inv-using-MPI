#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <unordered_set>
#include <mpi.h>

using namespace std;

// ------------------- Dense Generator -------------------
void generateDenseMatrix(vector<int>& mat, int n) {
    srand(time(NULL));
    for (int i = 0; i < n * n; ++i)
        mat[i] = rand() % 10;
}

// ------------------- Sparse Generator -------------------
void generateSparseMatrixAsDense(vector<int>& mat, int n, double sparsity) {
    srand(time(NULL));
    fill(mat.begin(), mat.end(), 0); // All zero first

    int non_zeros_per_row = static_cast<int>(n * sparsity);
    for (int i = 0; i < n; ++i) {
        unordered_set<int> cols;
        while (cols.size() < non_zeros_per_row) {
            int col = rand() % n;
            mat[i * n + col] = rand() % 10;
            cols.insert(col);
        }
    }
}

// ------------------- Multiply Block -------------------
vector<int> multiplyBlock(const vector<int>& A, const vector<int>& B, int n, int startRow, int endRow) {
    vector<int> local_result((endRow - startRow) * n, 0);

    for (int i = startRow; i < endRow; ++i) {
        for (int j = 0; j < n; ++j) {
            int sum = 0;
            for (int k = 0; k < n; ++k) {
                sum += A[i * n + k] * B[k * n + j];
            }
            local_result[(i - startRow) * n + j] = sum;
        }
    }
    return local_result;
}

// ------------------- Main -------------------
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n;
    char type;
    double sparsity = 0.0;

    if (rank == 0) {
        cout << "📦 Parallel Matrix Multiplication Program (MPI)\n";
        cout << "Choose matrix type:\n[d] Dense\n[s] Sparse\n> ";
        cin >> type;

        cout << "Enter size of square matrix (e.g., 3000): ";
        cin >> n;

        if (type == 's' || type == 'S') {
            cout << "Enter sparsity (0.0 to 1.0): ";
            cin >> sparsity;
        }
    }

    // Broadcast metadata
    MPI_Bcast(&type, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sparsity, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    vector<int> A(n * n, 0), B(n * n, 0);

    if (rank == 0) {
        if (type == 'd' || type == 'D') {
            generateDenseMatrix(A, n);
            generateDenseMatrix(B, n);
        } else {
            generateSparseMatrixAsDense(A, n, sparsity);
            generateSparseMatrixAsDense(B, n, sparsity);
        }
    }

    // One-shot broadcast of full matrices
    MPI_Bcast(A.data(), n * n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(B.data(), n * n, MPI_INT, 0, MPI_COMM_WORLD);

    // Work split
    int rows_per_proc = n / size;
    int start_row = rank * rows_per_proc;
    int end_row = (rank == size - 1) ? n : start_row + rows_per_proc;

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    // Local multiplication
    vector<int> local_result = multiplyBlock(A, B, n, start_row, end_row);

    double computation_end = MPI_Wtime();

    // Gather results
    vector<int> full_result;
    if (rank == 0)
        full_result.resize(n * n);

    MPI_Gather(local_result.data(), (end_row - start_row) * n, MPI_INT,
               full_result.data(), (end_row - start_row) * n, MPI_INT,
               0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();

    if (rank == 0) {
        double total_time = end_time - start_time;
        double comp_time = computation_end - start_time;
        double comm_time = total_time - comp_time;
        double comm_percent = (comm_time / total_time) * 100.0;

        if (type == 'd' || type == 'D')
            cout << "✅ Dense matrix multiplication completed.\n";
        else
            cout << "✅ Sparse matrix multiplication (via dense) completed.\n";

        cout << "⏱ Total Time: " << total_time << " seconds\n";
        cout << "🧠 Computation Time: " << comp_time << " seconds\n";
        cout << "📡 Communication Time: " << comm_time << " seconds\n";
        cout << "📊 Communication %: " << comm_percent << "%\n";
    }

    MPI_Finalize();
    return 0;
}
