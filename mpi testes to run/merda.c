#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define DIM_X 6
#define DIM_Y 6
#define DIM_Z 6

char*** create3DArray(int dim_x, int dim_y, int dim_z) {
    char*** array = (char***)malloc(dim_x * sizeof(char**));
    for (int i = 0; i < dim_x; ++i) {
        array[i] = (char**)malloc(dim_y * sizeof(char*));
        for (int j = 0; j < dim_y; ++j) {
            array[i][j] = (char*)malloc(dim_z * sizeof(char));
        }
    }
    return array;
}

void fillArray(char*** array, int dim_x, int dim_y, int dim_z) {
    for (int i = 0; i < dim_x; ++i) {
        for (int j = 0; j < dim_y; ++j) {
            for (int k = 0; k < dim_z; ++k) {
                array[i][j][k] = 'A' + (i + j + k) % 26; // Filling array with characters
            }
        }
    }
}

void printLocalArray(char*** local_array, int local_dim_x, int local_dim_y, int local_dim_z) {
    for (int i = 0; i < local_dim_x; ++i) {
        for (int j = 0; j < local_dim_y; ++j) {
            for (int k = 0; k < local_dim_z; ++k) {
                printf("%c ", local_array[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int dims[3] = {0, 0, 0}; // To store dimensions of the cartesian grid
    int periods[3] = {0, 0, 0}; // To store whether the grid is periodic or not
    MPI_Dims_create(size, 3, dims);

    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart_comm);

    int coords[3];
    MPI_Cart_coords(cart_comm, rank, 3, coords);

    int local_dim_x = DIM_X / dims[0];
    int local_dim_y = DIM_Y / dims[1];
    int local_dim_z = DIM_Z / dims[2];

    char ***local_array = create3DArray(local_dim_x, local_dim_y, local_dim_z);

    if (rank == 0) {
        char ***global_array = create3DArray(DIM_X, DIM_Y, DIM_Z);
        fillArray(global_array, DIM_X, DIM_Y, DIM_Z);
        // Scatter the global array among processes
        MPI_Scatter(global_array[0][0], local_dim_x * local_dim_y * local_dim_z, MPI_CHAR,
                    local_array[0][0], local_dim_x * local_dim_y * local_dim_z, MPI_CHAR,
                    0, MPI_COMM_WORLD);
        // Printing local arrays
        printf("Process %d, Coordinates (%d, %d, %d):\n", rank, coords[0], coords[1], coords[2]);
        printLocalArray(local_array, local_dim_x, local_dim_y, local_dim_z);

        // Clean up
        for (int i = 0; i < DIM_X; ++i) {
            for (int j = 0; j < DIM_Y; ++j) {
                free(global_array[i][j]);
            }
            free(global_array[i]);
        }
        free(global_array);
    } else {
        // Scatter the global array among processes
        MPI_Scatter(NULL, local_dim_x * local_dim_y * local_dim_z, MPI_CHAR,
                    local_array[0][0], local_dim_x * local_dim_y * local_dim_z, MPI_CHAR,
                    0, MPI_COMM_WORLD);
        // Printing local arrays
        printf("Process %d, Coordinates (%d, %d, %d):\n", rank, coords[0], coords[1], coords[2]);
        printLocalArray(local_array, local_dim_x, local_dim_y, local_dim_z);
    }

    // Clean up
    for (int i = 0; i < local_dim_x; ++i) {
        for (int j = 0; j < local_dim_y; ++j) {
            free(local_array[i][j]);
        }
        free(local_array[i]);
    }
    free(local_array);

    MPI_Finalize();
    return 0;
}
