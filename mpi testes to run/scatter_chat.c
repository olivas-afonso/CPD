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

void reconstructMatrix(int *flatMatrix, int received_matrix[X_DIM][Y_DIM]) {
    for (int i = 0; i < X_DIM; i++) {
        for (int j = 0; j < Y_DIM; j++) {
            received_matrix[i][j] = flatMatrix[i * Y_DIM + j];
        }
    }
}

void printMatrix(int received_matrix[X_DIM][Y_DIM], int rank) {
    for (int i = 0; i < X_DIM; i++) {
        for (int j = 0; j < Y_DIM; j++) {
            printf("Process %d: matrix[%d][%d] = %d\n", rank, i, j, received_matrix[i][j]);
        }
    }
    printf("\n");
}

int main(int argc, char **argv) {
    int rank, size;
    int sendcounts[X_DIM], displs[X_DIM];
    int flatMatrix[X_DIM * Y_DIM];
    int recvBuffer[X_DIM * Y_DIM]; // Separate receive buffer for each process
    int received_matrix[X_DIM][Y_DIM]; // Matrix to hold received data

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Datatype matrix3d_type;
    MPI_Type_contiguous(X_DIM * Y_DIM * Z_DIM, MPI_INT, &matrix3d_type);
    MPI_Type_commit(&matrix3d_type);

    // Prepare data to scatter
    if (rank == 0) {
        // Initialize the matrix
        Matrix3D matrix;
        for (int i = 0; i < X_DIM; i++) {
            for (int j = 0; j < Y_DIM; j++) {
                for (int k = 0; k < Z_DIM; k++) {
                    matrix.data[i][j][k] = (i * Y_DIM * Z_DIM) + (j * Z_DIM) + k; // Ensure values are within manageable range
                }
            }
        }

        // Print the initial matrix
        printf("Initial matrix (rank 0):\n");
        for (int i = 0; i < X_DIM; i++) {
            for (int j = 0; j < Y_DIM; j++) {
                printf("%d ", matrix.data[i][j][0]);
            }
            printf("\n");
        }
        printf("\n");

        // Flatten the matrix
        flattenMatrix(&matrix, flatMatrix);

        // Prepare sendcounts and displs for scatterv
        for (int i = 0; i < X_DIM; i++) {
            sendcounts[i] = Y_DIM;
            displs[i] = i * Y_DIM;
        }
        
        // Print sendcounts and displs for debugging
        printf("Send counts (rank 0): ");
        for (int i = 0; i < X_DIM; i++) {
            printf("%d ", sendcounts[i]);
        }
        printf("\n");
        printf("Displacements (rank 0): ");
        for (int i = 0; i < X_DIM; i++) {
            printf("%d ", displs[i]);
        }
        printf("\n");
    }

    // Ensure all processes reach this point before proceeding
    MPI_Barrier(MPI_COMM_WORLD);

    // Scatter the matrix to all processes
    MPI_Scatterv(flatMatrix, sendcounts, displs, MPI_INT, recvBuffer, X_DIM * Y_DIM, MPI_INT, 0, MPI_COMM_WORLD);

    // Print received data for debugging
    MPI_Gather(recvBuffer, X_DIM * Y_DIM, MPI_INT, flatMatrix, X_DIM * Y_DIM, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("Received matrices:\n");
        for (int p = 0; p < size; p++) {
            reconstructMatrix(&flatMatrix[p * X_DIM * Y_DIM], received_matrix);
            printf("Received matrix (rank %d):\n", p);
            printMatrix(received_matrix, p);
        }
    }

    MPI_Type_free(&matrix3d_type);
    MPI_Finalize();
    return 0;
}
