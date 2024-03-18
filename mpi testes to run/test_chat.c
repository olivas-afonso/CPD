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

    // Determine layer number for this process
    int layer_num = coords[2];

    // Create a 3D array to hold the layer of the grid for each process
    int layer[NX][NY][NZ];
    
    // Fill out the layer of the grid with example values
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                layer[i][j][k] = layer_num + 1; // Fill with the layer number
            }
        }
    }

    // Print out the layer of the grid for each process
    printf("Rank %d: Layer %d\n", rank, layer_num);
    printLayer(layer, layer_num);

    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
    return 0;
}
