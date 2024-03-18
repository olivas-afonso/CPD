#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define X_DIM 3
#define Y_DIM 3
#define Z_DIM 3

typedef struct {
    int data[X_DIM][Y_DIM][Z_DIM];
} Matrix3D;

int main(int argc, char **argv) {
    int rank, size;
    int sendcounts, displs;
    Matrix3D matrix, received_matrix;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Datatype matrix3d_type;
    MPI_Type_contiguous(X_DIM * Y_DIM * Z_DIM, MPI_INT, &matrix3d_type);
    MPI_Type_commit(&matrix3d_type);

    // Prepare data to scatter
    if (rank == 0) {
        // Initialize the matrix
        for (int i = 0; i < X_DIM; i++) {
            for (int j = 0; j < Y_DIM; j++) {
                for (int k = 0; k < Z_DIM; k++) {
                    matrix.data[i][j][k] = i + j + k;
                }
            }
        }

        // Prepare sendcounts and displs for scatterv
        sendcounts = X_DIM * Y_DIM * Z_DIM;
        displs = 0;
    }

    // Scatter the matrix to all processes
    MPI_Scatterv(&matrix, &sendcounts, &displs, matrix3d_type, &received_matrix, 1, matrix3d_type, 0, MPI_COMM_WORLD);

    // Each process prints the received data
    for (int x = 0; x < X_DIM; x++) {
        for (int y = 0; y < Y_DIM; y++) {
            for (int z = 0; z < Z_DIM; z++) {
                printf("Process %d: received matrix[%d][%d][%d] = %d\n", rank, x, y, z, received_matrix.data[x][y][z]);
            }
        }
    }

    MPI_Type_free(&matrix3d_type);
    MPI_Finalize();
    return 0;
}
