#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define NX 12
#define NY 8
#define NZ 6

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

    // Each process gets responsibility for one layer
    int local_size_z = NZ / dims[2];
    int start_z = coords[2] * local_size_z;
    int end_z = (coords[2] + 1) * local_size_z;

    // Fill out the grid with example values
    int grid[NX][NY][NZ];
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                grid[i][j][k] = (i * NY * NZ + j * NZ + k) + 1; // Unique values for each cell
            }
        }
    }

    // Print out the filled grid
    printf("Rank %d: Z: %d-%d\n", rank, start_z, end_z);

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = start_z; k < end_z; k++) {
                printf("Rank %d received value %d at [%d][%d][%d]\n", rank, grid[i][j][k], i, j, k);
            }
        }
    }

    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
    return 0;
}