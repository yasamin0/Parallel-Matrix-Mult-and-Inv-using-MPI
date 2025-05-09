#include <iostream>     // For input/output
#include <vector>       // For using std::vector to represent matrices
#include <cstdlib>      // For rand(), srand()
#include <ctime>        // For time-based random seed
#include <mpi.h>        // For MPI parallel programming

using namespace std;

// ------------------ Gauss-Jordan Row Operation ------------------
// Performs distributed Gauss-Jordan elimination for matrix inversion
void gaussJordan(vector<vector<double>>& matrix, vector<vector<double>>& identity, int n, int rank, int size) {
    for (int i = 0; i < n; i++) {
        int owner = i * size / n;  // Determine which rank owns the i-th row

        // The owner of the pivot row performs normalization
        if (rank == owner) {
            double pivot = matrix[i][i];
            if (pivot == 0.0) {
                if (rank == 0) cerr << "Matrix is singular or nearly singular.\n";
                MPI_Abort(MPI_COMM_WORLD, 1);  // Abort if matrix is not invertible
            }
            // Normalize pivot row
            for (int j = 0; j < n; j++) {
                matrix[i][j] /= pivot;
                identity[i][j] /= pivot;
            }
        }

        // Broadcast the normalized pivot row to all processes
        MPI_Bcast(matrix[i].data(), n, MPI_DOUBLE, owner, MPI_COMM_WORLD);
        MPI_Bcast(identity[i].data(), n, MPI_DOUBLE, owner, MPI_COMM_WORLD);

        // Eliminate pivot column from other rows (that this rank owns)
        for (int k = 0; k < n; ++k) {
            if (k == i || (k * size / n) != rank) continue;  // Skip the pivot row and rows not owned by this rank
            // akhe k hamun satr mehvari bashe ya age in satr be pardazande feli rabti nadare nabayad bahash kari dashte bashim 
            double factor = matrix[k][i];
            for (int j = 0; j < n; ++j) {
                matrix[k][j] -= factor * matrix[i][j];
                identity[k][j] -= factor * identity[i][j];
            }
        }
    }
}

// ------------------ Main Program ------------------
int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);  // Initialize MPI environment

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the current rank (process ID)
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the total number of processes

    int n;  // Size of the matrix
    vector<vector<double>> matrix, identity;

    // Only rank 0 initializes the matrix
    if (rank == 0) {
        cout << "Parallel Matrix Inversion (Gauss-Jordan, MPI)\n";
        cout << "Enter size of square matrix (e.g., 1000): ";
        cin >> n;

        srand(time(NULL)); //تابع تولید عدد تصادفی رو بر اساس زمان فعلی تنظیم کن
        matrix.resize(n, vector<double>(n));
        identity.resize(n, vector<double>(n, 0.0));

        // Initialize matrix with random values and identity matrix
        for (int i = 0; i < n; ++i) {
            identity[i][i] = 1.0;
            for (int j = 0; j < n; ++j)
                matrix[i][j] = rand() % 10 + 1;  // Avoid zeros on diagonal
        }
    }

    // Broadcast matrix size to all processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Other ranks misazan memory for matrices and initialize identity
    if (rank != 0) {
        matrix.resize(n, vector<double>(n));
        identity.resize(n, vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i)
            identity[i][i] = 1.0;
    }

    // Broadcast matrix data row by row from rank 0 to all processes
    for (int i = 0; i < n; ++i) {
        MPI_Bcast(matrix[i].data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(identity[i].data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Synchronize before starting the timing
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();  // Start total timing

    // Perform parallel Gauss-Jordan elimination
    gaussJordan(matrix, identity, n, rank, size);

    double computation_end = MPI_Wtime();  // End of computation

    MPI_Barrier(MPI_COMM_WORLD);  // Ensure all processes finish
    double end_time = MPI_Wtime();  // End of total timing

    // Only rank 0 reports performance metrics
    if (rank == 0) {
        double total_time = end_time - start_time;
        double comp_time = computation_end - start_time;
        double comm_time = total_time - comp_time;
        double comm_percent = (comm_time / total_time) * 100.0;

        cout << "Matrix inversion completed.\n";
        cout << "Total Time: " << total_time << " seconds\n";
        cout << "Computation Time: " << comp_time << " seconds\n";
        cout << "Communication Time: " << comm_time << " seconds\n";
        cout << "Communication %: " << comm_percent << "%\n";
    }

    MPI_Finalize();  // Finalize MPI
    return 0;
}
