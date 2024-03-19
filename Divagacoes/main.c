#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define WIDTH 12
#define HEIGHT 12
#define DEPTH 12

// Function to initialize the 3D matrix
void initializeMatrix(int matrix[WIDTH][HEIGHT][DEPTH]) {
    for (int x = 0; x < WIDTH; ++x) {
        for (int y = 0; y < HEIGHT; ++y) {
            for (int z = 0; z < DEPTH; ++z) {
                matrix[x][y][z] = x + y + z; // Example initialization
            }
        }
    }
}

// Function to print the 3D matrix
void printMatrix(int matrix[WIDTH][HEIGHT][DEPTH]) {
    for (int z = 0; z < DEPTH; ++z) {
        printf("Layer %d:\n", z);
        for (int x = 0; x < WIDTH; ++x) {
            for (int y = 0; y < HEIGHT; ++y) {
                printf("%d ", matrix[x][y][z]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

int main(int argc, char** argv) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Calculate the size of each sub-matrix
    int subWidth = WIDTH / size;
    int subHeight = HEIGHT / size;
    int subDepth = DEPTH / size;

    // Check if the number of processes is compatible with the matrix dimensions
    if (WIDTH % size != 0 || HEIGHT % size != 0 || DEPTH % size != 0) {
        if (rank == 0) {
            printf("Error: Number of processes must evenly divide matrix dimensions.\n");
        }
        MPI_Finalize();
        return 1;
    }

    int subMatrix[subWidth][subHeight][subDepth];

    // Initialize the 3D matrix on process 0
    if (rank == 0) {
        int matrix[WIDTH][HEIGHT][DEPTH];
        initializeMatrix(matrix);
        printf("Initial matrix:\n");
        printMatrix(matrix);

        // Scatter the matrix data to all processes
        MPI_Scatter(matrix, subWidth * subHeight * subDepth, MPI_INT, subMatrix,
                    subWidth * subHeight * subDepth, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        // Receive the sub-matrix data from process 0
        MPI_Scatter(NULL, subWidth * subHeight * subDepth, MPI_INT, subMatrix,
                    subWidth * subHeight * subDepth, MPI_INT, 0, MPI_COMM_WORLD);
    }

    // Each process prints its own sub-matrix
    printf("Rank %d sub-matrix:\n", rank);
    printMatrix(subMatrix);

    MPI_Finalize();
    return 0;
}
