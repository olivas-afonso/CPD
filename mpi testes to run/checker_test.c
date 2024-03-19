#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N_X 100
#define N_Y 100
#define N_Z 100

void printSubgrid(char ***subgrid, int subgrid_size_x, int subgrid_size_y, int subgrid_size_z) {
    for (int i = 0; i < subgrid_size_x; i++) {
        for (int j = 0; j < subgrid_size_y; j++) {
            for (int k = 0; k < subgrid_size_z; k++) {
                printf("%c ", subgrid[i][j][k]);
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

    // Create subgrid for this process only if it has a valid portion of the grid
    char ***subgrid = NULL;
    if (subgrid_size_x > 0 && subgrid_size_y > 0 && subgrid_size_z > 0) {
        subgrid = (char ***)malloc(subgrid_size_x * sizeof(char **));
        for (int i = 0; i < subgrid_size_x; i++) {
            subgrid[i] = (char **)malloc(subgrid_size_y * sizeof(char *));
            for (int j = 0; j < subgrid_size_y; j++) {
                subgrid[i][j] = (char *)malloc(subgrid_size_z * sizeof(char));
            }
        }

        // Fill subgrid with some data (for demonstration)
        for (int i = 0; i < subgrid_size_x; i++) {
            for (int j = 0; j < subgrid_size_y; j++) {
                for (int k = 0; k < subgrid_size_z; k++) {
                    subgrid[i][j][k] = 'A' + rank; // Filling with characters starting from 'A' + rank
                }
            }
        }
    }

    // Allocate memory for gathered_subgrids for all processes
    char ***gathered_subgrids = (char ***)malloc(N_X * sizeof(char **));
    for (int i = 0; i < N_X; i++) {
        gathered_subgrids[i] = (char **)malloc(N_Y * sizeof(char *));
        for (int j = 0; j < N_Y; j++) {
            gathered_subgrids[i][j] = (char *)malloc(N_Z * sizeof(char));
        }
    }

    // Gather subgrids from all processes to process 0
    int gather_status = MPI_Gather((subgrid != NULL) ? &(subgrid[0][0][0]) : NULL, 
                                   (subgrid != NULL) ? (subgrid_size_x * subgrid_size_y * subgrid_size_z) : 0, MPI_CHAR, 
                                   &(gathered_subgrids[0][0][0]), (subgrid != NULL) ? (subgrid_size_x * subgrid_size_y * subgrid_size_z) : 0, MPI_CHAR, 0, MPI_COMM_WORLD);

    // Process 0 prints the gathered subgrid
    if (rank == 0) {
        printf("Gathered Subgrid:\n");
        printSubgrid(gathered_subgrids, N_X, N_Y, N_Z);
    }

    // Free memory
    if (subgrid != NULL) {
        for (int i = 0; i < subgrid_size_x; i++) {
            for (int j = 0; j < subgrid_size_y; j++) {
                free(subgrid[i][j]);
            }
            free(subgrid[i]);
        }
        free(subgrid);
    }

    for (int i = 0; i < N_X; i++) {
        for (int j = 0; j < N_Y; j++) {
            free(gathered_subgrids[i][j]);
        }
        free(gathered_subgrids[i]);
    }
    free(gathered_subgrids);

    MPI_Finalize();

    // Print diagnostic information
    if (gather_status != MPI_SUCCESS) {
        fprintf(stderr, "Error: MPI_Gather failed with status %d\n", gather_status);
        return 1;
    }

    return 0;
}
