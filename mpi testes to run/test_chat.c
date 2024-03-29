#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define NX 6
#define NY 6
#define NZ 6

void printLayer(int layer[NX][NY][NZ], int layer_num, int rank) {
    printf("Rank %d: Layer %d\n", rank, layer_num);
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
    int periods[3] = {1, 1, 1};
    int coords[3];
    MPI_Comm cart_comm;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    printf("SIZE: %d\n", size);
    // Create 3D Cartesian communicator
    MPI_Dims_create(size, 3, dims);
    
    printf("DIMS %d %d %d:\n", dims[0], dims[1], dims[2]);
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &cart_comm);
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

    // Synchronize the output
    MPI_Barrier(MPI_COMM_WORLD);

    // Only let process 0 print the initial grid
    if (rank == 0) {
        printf("Initial Grid:\n");
        for (int k = 0; k < NZ; k++) {
            printLayer(layer, k, rank);
        }
    }

    // Synchronize the output
    MPI_Barrier(MPI_COMM_WORLD);

    // Print out the layer of the grid for each process one at a time
    for (int i = 0; i < size; i++) {
        if (rank == i) {
            printf("Rank %d is printing...\n", rank);
            for (int k = start_layer; k < end_layer; k++) {
                MPI_Barrier(MPI_COMM_WORLD); // Synchronize before receiving
                MPI_Bcast(&layer, NX * NY * NZ, MPI_INT, i, MPI_COMM_WORLD); // Broadcast the layer
                printLayer(layer, k, rank); // Print the received layer
                if (k != end_layer - 1) {
                    printf("\n"); // Print new line unless it's the last layer
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
    return 0;
}
