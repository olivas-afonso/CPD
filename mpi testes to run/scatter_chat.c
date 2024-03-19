#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define X_DIM 3
#define Y_DIM 3



unsigned int seed;

typedef struct {
    int data***;
} Matrix3D;

void flattenMatrix(Matrix3D *matrix, int *flatMatrix, int Z_DIM) {
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

    int number_of_gens,  number_of_cells;
    float density;

    number_of_gens = atoi (argv[1]);
    number_of_cells = atoi (argv[2]);
    density = atof (argv[3]);
    seed = atoi (argv[4]);

    int Z_DIM = number_of_cells;
    
    int rank, size;
    int sendcounts[X_DIM], displs[X_DIM];
    int flatMatrix[X_DIM * Y_DIM * Z_DIM]={0};
    int recvBuffer[X_DIM * Y_DIM]={0}; // Separate receive buffer for each process
    int received_matrix[X_DIM][Y_DIM];
    Matrix3D matrix;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Datatype matrix3d_type;
    MPI_Type_contiguous(X_DIM * Y_DIM * Z_DIM, MPI_INT, &matrix3d_type);
    MPI_Type_commit(&matrix3d_type);

    

    for (int i = 0; i < X_DIM; i++) {
        for (int j = 0; j < Y_DIM; j++) {

                received_matrix[i][j] = 0;
   
        }
    }

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
            for(int aux_x4=0; aux_x4<X_DIM; aux_x4++)
            {
                for(int aux_y4=0; aux_y4<Y_DIM; aux_y4++)
                {
                    for(int aux_z4=0; aux_z4<Z_DIM; aux_z4++)
                    {
                        printf("%d ", matrix.data[aux_x4][aux_y4][aux_z4]);
                    }
                    printf("\n");
                }
                printf("\n");
            }

        // Flatten the matrix
        flattenMatrix(&matrix, flatMatrix, Z_DIM);

        // Prepare sendcounts and displs for scatterv
        for (int i = 0; i < X_DIM; i++) {
            sendcounts[i] = Y_DIM * Z_DIM;
            displs[i] = i * Y_DIM * Z_DIM;
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
    MPI_Scatterv(flatMatrix, sendcounts, displs, MPI_INT, recvBuffer, Y_DIM * Z_DIM, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    // Print received data for debugging
    printf("Received data (rank %d): ", rank);
    for (int i = 0; i < Y_DIM * Z_DIM; i++) {
        printf("%d ", recvBuffer[i]);
    }
    printf("\n");

    // Ensure all processes receive their data before proceeding
    MPI_Barrier(MPI_COMM_WORLD);

    // Reconstruct the matrix
    reconstructMatrix(recvBuffer, received_matrix);

    MPI_Barrier(MPI_COMM_WORLD);
    // Print the received matrix for each process
    printf("Received matrix (rank %d):\n", rank);
    printMatrix(received_matrix, rank);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Type_free(&matrix3d_type);
    MPI_Finalize();
    return 0;
}
