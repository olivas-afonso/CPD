#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N_X 100
#define N_Y 100
#define N_Z 100

void printGrid(char ***grid, int dim_x, int dim_y, int dim_z) {
    for (int i = 0; i < dim_x; i++) {
        for (int j = 0; j < dim_y; j++) {
            for (int k = 0; k < dim_z; k++) {
                printf("%c ", grid[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Calculate subgrid size for each dimension
    int subgrid_size_x = N_X / size; // Divide N_X by the number of processes
    int subgrid_size_y = N_Y / size; // Divide N_Y by the number of processes
    int subgrid_size_z = N_Z / size; // Divide N_Z by the number of processes

    // Calculate start and end indices along the X-axis for each process
    int start_x = rank * subgrid_size_x;
    int end_x = start_x + subgrid_size_x;

    // Create subgrid for each process
    char ***subgrid = (char ***)malloc(subgrid_size_x * sizeof(char **));
    for (int i = 0; i < subgrid_size_x; i++) {
        subgrid[i] = (char **)malloc(subgrid_size_y * sizeof(char *));
        for (int j = 0; j < subgrid_size_y; j++) {
            subgrid[i][j] = (char *)malloc(subgrid_size_z * sizeof(char));
        }
    }

    // Perform computation on the subgrid (for demonstration, let's fill it with process rank)
    for (int i = 0; i < subgrid_size_x; i++) {
        for (int j = 0; j < subgrid_size_y; j++) {
            for (int k = 0; k < subgrid_size_z; k++) {
                subgrid[i][j][k] = 'A' + rank; // Filling with characters starting from 'A' + rank
            }
        }
    }

    // Gather results from all processes
    char ***gathered_subgrids = NULL;
    if (rank == 0) {
        gathered_subgrids = (char ***)malloc(N_X * sizeof(char **));
        for (int i = 0; i < N_X; i++) {
            gathered_subgrids[i] = (char **)malloc(N_Y * sizeof(char *));
            for (int j = 0; j < N_Y; j++) {
                gathered_subgrids[i][j] = (char *)malloc(N_Z * sizeof(char));
            }
        }
    }
    printf("TESTING \n");
    MPI_Gather(&(subgrid[0][0][0]), subgrid_size_x * subgrid_size_y * subgrid_size_z, MPI_CHAR, 
               &(gathered_subgrids[0][0][0]), subgrid_size_x * subgrid_size_y * subgrid_size_z, MPI_CHAR, 0, MPI_COMM_WORLD);

    // Process 0 aggregates all subgrids
    if (rank == 0) {
        // Print the gathered subgrid
        printf("Gathered Subgrid:\n");
        printGrid(gathered_subgrids, N_X, N_Y, N_Z);
    }

    // Free memory
    for (int i = 0; i < subgrid_size_x; i++) {
        for (int j = 0; j < subgrid_size_y; j++) {
            free(subgrid[i][j]);
        }
        free(subgrid[i]);
    }
    free(subgrid);

    if (rank == 0) {
        for (int i = 0; i < N_X; i++) {
            for (int j = 0; j < N_Y; j++) {
                free(gathered_subgrids[i][j]);
            }
            free(gathered_subgrids[i]);
        }
        free(gathered_subgrids);
    }

    MPI_Finalize();
    return 0;
}
 