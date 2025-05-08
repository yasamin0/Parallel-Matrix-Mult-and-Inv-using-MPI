#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <mpi.h>

using namespace std;

// ------------------ Gauss-Jordan Row Operation ------------------
void gaussJordan(vector<vector<double>>& matrix, vector<vector<double>>& identity, int n, int rank, int size) {
    for (int i = 0; i < n; i++) {
        int owner = i * size / n;
        if (rank == owner) {
            double pivot = matrix[i][i];
            if (pivot == 0.0) {
                if (rank == 0) cerr << "Matrix is singular or nearly singular.\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            for (int j = 0; j < n; j++) {
                matrix[i][j] /= pivot;
                identity[i][j] /= pivot;
            }
        }

        MPI_Bcast(matrix[i].data(), n, MPI_DOUBLE, owner, MPI_COMM_WORLD);
        MPI_Bcast(identity[i].data(), n, MPI_DOUBLE, owner, MPI_COMM_WORLD);

        for (int k = 0; k < n; ++k) {
            if (k == i || (k * size / n) != rank) continue;

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
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n;
    vector<vector<double>> matrix, identity;

    if (rank == 0) {
        cout << "🧮 Parallel Matrix Inversion (Gauss-Jordan, MPI)\n";
        cout << "Enter size of square matrix (e.g., 1000): ";
        cin >> n;

        srand(time(NULL));
        matrix.resize(n, vector<double>(n));
        identity.resize(n, vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i) {
            identity[i][i] = 1.0;
            for (int j = 0; j < n; ++j)
                matrix[i][j] = rand() % 10 + 1;
        }
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        matrix.resize(n, vector<double>(n));
        identity.resize(n, vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i)
            identity[i][i] = 1.0;
    }

    for (int i = 0; i < n; ++i) {
        MPI_Bcast(matrix[i].data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(identity[i].data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    gaussJordan(matrix, identity, n, rank, size);

    double computation_end = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();

    if (rank == 0) {
        double total_time = end_time - start_time;
        double comp_time = computation_end - start_time;
        double comm_time = total_time - comp_time;
        double comm_percent = (comm_time / total_time) * 100.0;

        cout << "✅ Matrix inversion completed.\n";
        cout << "⏱ Total Time: " << total_time << " seconds\n";
        cout << "🧠 Computation Time: " << comp_time << " seconds\n";
        cout << "📡 Communication Time: " << comm_time << " seconds\n";
        cout << "📊 Communication %: " << comm_percent << "%\n";
    }

    MPI_Finalize();
    return 0;
}
