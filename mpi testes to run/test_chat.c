#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define NX 4
#define NY 4
#define NZ 4

void printLayer(int layer[NX][NY][NZ], int layer_num) {
    printf("Layer %d: \n", layer_num);
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            printf("%2d ", layer[i][j][layer_num]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char **argv) {
    int rank, size;
    int dims[3] = {0, 0, 0};
    int periods[3] = {0, 0, 0};
    int coords[3];
    MPI_Comm cart_comm;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Create 3D Cartesian communicator
    MPI_Dims_create(size, 3, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 3, coords);

    // Each process gets responsibility for a layer of the grid
    int layers_per_process = NZ / dims[2];
    int remainder = NZ % dims[2]; // Check for any remaining layers
    int start_layer = coords[2] * layers_per_process + (coords[2] < remainder ? coords[2] : remainder);
    int end_layer = start_layer + layers_per_process + (coords[2] < remainder ? 1 : 0);
    
    // Create a 3D array to hold the layer of the grid for each process
    int layer[NX][NY][NZ];
    
    // Fill out the layer of the grid with example values
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                layer[i][j][k] = k + 1; // Fill with the number of the layer
            }
        }
    }

    // Print out the initial grid (only by process 0)
    if (rank == 0) {
        printf("Initial Grid:\n");
        for (int k = 0; k < NZ; k++) {
            printLayer(layer, k);
        }
    }

    // Wait for process 0 to finish printing initial grid
    MPI_Barrier(MPI_COMM_WORLD);

    // Print out the layer of the grid for each process
    printf("Rank %d: Layer %d-%d\n", rank, start_layer, end_layer - 1);
    for (int k = start_layer; k < end_layer; k++) {
        printLayer(layer, k);
        if (k != end_layer - 1) {
            printf("\n"); // Print new line unless it's the last layer
        }
    }

    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
    return 0;
}
