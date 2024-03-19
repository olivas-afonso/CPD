#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define X_DIM 3
#define Y_DIM 3
#define Z_DIM 3

typedef struct {
    int data[X_DIM][Y_DIM][Z_DIM];
} Matrix3D;

void flattenMatrix(Matrix3D *matrix, int *flatMatrix) {
    for (int i = 0; i < X_DIM; i++) {
        for (int j = 0; j < Y_DIM; j++) {
            for (int k = 0; k < Z_DIM; k++) {
                flatMatrix[i * Y_DIM * Z_DIM + j * Z_DIM + k] = matrix->data[i][j][k];
            }
        }
    }
}

void reconstructMatrix(int *flatMatrix, Matrix3D *matrix) {
    for (int i = 0; i < X_DIM; i++) {
        for (int j = 0; j < Y_DIM; j++) {
            for (int k = 0; k < Z_DIM; k++) {
                matrix->data[i][j][k] = flatMatrix[i * Y_DIM * Z_DIM + j * Z_DIM + k];
            }
        }
    }
}

void printMatrix(Matrix3D *matrix, int rank) {
    for (int i = 0; i < X_DIM; i++) {
        for (int j = 0; j < Y_DIM; j++) {
            for (int k = 0; k < Z_DIM; k++) {
                printf("Process %d: matrix[%d][%d][%d] = %d\n", rank, i, j, k, matrix->data[i][j][k]);
            }
        }
    }
    printf("\n");
}

int main(int argc, char **argv) {
    int rank, size;
    int sendcounts[X_DIM], displs[X_DIM];
    int flatMatrix[X_DIM * Y_DIM * Z_DIM];
    int recvBuffer[Y_DIM * Z_DIM]; // Separate receive buffer for each process
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
                    matrix.data[i][j][k] = (i * Y_DIM * Z_DIM) + (j * Z_DIM) + k; // Ensure values are within manageable range
                }
            }
        }

        // Print the initial matrix
        printf("Initial matrix (rank 0):\n");
        printMatrix(&matrix, rank);

        // Flatten the matrix
        flattenMatrix(&matrix, flatMatrix);

        // Prepare sendcounts and displs for scatterv
        for (int i = 0; i < X_DIM; i++) {
            sendcounts[i] = Y_DIM * Z_DIM;
            displs[i] = i * Y_DIM * Z_DIM;
        }
    }

    // Ensure all processes reach this point before proceeding
    MPI_Barrier(MPI_COMM_WORLD);

    // Scatter the matrix to all processes
    MPI_Scatterv(flatMatrix, sendcounts, displs, MPI_INT, recvBuffer, Y_DIM * Z_DIM, MPI_INT, 0, MPI_COMM_WORLD);

    // Ensure all processes receive their data before proceeding
    MPI_Barrier(MPI_COMM_WORLD);

    // Reconstruct the matrix
    reconstructMatrix(recvBuffer, &received_matrix);

    // Print the received matrix for each process
    printf("Received matrix (rank %d):\n", rank);
    printMatrix(&received_matrix, rank);

    MPI_Type_free(&matrix3d_type);
    MPI_Finalize();
    return 0;
}
