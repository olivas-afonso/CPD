#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define NX 12
#define NY 8
#define NZ 6

void printLayer(int layer[NX][NY][NZ], int start_layer, int end_layer) {
    for (int k = start_layer; k < end_layer; k++) {
        printf("Layer %d: \n", k);
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                printf("%2d ", layer[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {
    int rank, size;
    int dims[3] = {0, 0, 0};
    int periods[3] = {0, 0, 0};
    int coords[3];
    MPI_Comm cart_comm;
    int remain_dims[3] = {1, 1, 0}; // Do not wrap around in the Z dimension

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Create 3D Cartesian communicator
    MPI_Dims_create(size, 3, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 3, coords);

    // Each process gets responsibility for a layer of the grid
    int layers_per_process = NZ / dims[2];
    int start_layer = coords[2] * layers_per_process;
    int end_layer = (coords[2] + 1) * layers_per_process;

    // Create a 3D array to hold the layer of the grid for each process
    int layer[NX][NY][NZ];
    
    // Fill out the layer of the grid with example values
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = start_layer; k < end_layer; k++) {
                layer[i][j][k] = k + 1; // Fill with the number of the layer
            }
        }
    }

    // Print out the initial grid
    if (rank == 0) {
        printf("Initial Grid:\n");
        for (int k = 0; k < NZ; k++) {
            printf("Layer %d:\n", k);
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    printf("%2d ", k + 1);
                }
                printf("\n");
            }
            printf("\n");
        }
    }

    // Print out the layer of the grid for each process
    printf("Rank %d: Layer %d-%d\n", rank, start_layer, end_layer);
    printLayer(layer, start_layer, end_layer);

    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
    return 0;
}
