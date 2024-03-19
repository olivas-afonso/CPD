#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define WIDTH 12
#define HEIGHT 12
#define DEPTH 12

int valor = 0;

// Function to initialize the 3D matrix
void initializeMatrix(int matrix[WIDTH][HEIGHT][DEPTH]) {
    for (int x = 0; x < WIDTH; ++x) {
        for (int y = 0; y < HEIGHT; ++y) {
            for (int z = 0; z < DEPTH; ++z) {
                matrix[x][y][z] = valor; // Example initialization
                ++valor;
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

    // Calculate the remaining elements if dimensions are not evenly divisible
    int remainderWidth = WIDTH % size;
    int remainderHeight = HEIGHT % size;
    int remainderDepth = DEPTH % size;

    // Adjust the sub-matrix size if there's a remainder
    if (rank < remainderWidth) {
        subWidth++;
    }
    if (rank >= size - remainderWidth) {
        subWidth--;
    }

    if (rank < remainderHeight) {
        subHeight++;
    }
    if (rank >= size - remainderHeight) {
        subHeight--;
    }

    if (rank < remainderDepth) {
        subDepth++;
    }
    if (rank >= size - remainderDepth) {
        subDepth--;
    }

    // Allocate memory for sub-matrix
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
