#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define X_DIM 3
#define Y_DIM 3
#define N_SPECIES 9



unsigned int seed;

typedef struct {
    char *** data;
} Matrix3D;

Matrix3D matrix;

void init_r4uni(int input_seed)
{
    seed = input_seed + 987654321;
}

float r4_uni()
{
    int seed_in = seed;

    seed ^= (seed << 13);
    seed ^= (seed >> 17);
    seed ^= (seed << 5);

    return 0.5 + 0.2328306e-09 * (seed_in + (int) seed);
}

char *** allocated_matrix(int N, float density, int input_seed)
{

}


char ***gen_initial_grid(int N, float density, int input_seed)
{
    int x, y, z;
    
    init_r4uni(input_seed);
    
//alocacao da memeoria dinamica, alocando primeiro um apontador triplo o que corresponde a uma dimensao do cubo  ~
    matrix.data = (char ***) malloc(N * sizeof(char **));
    if(matrix.data == NULL) {
        printf("Failed to allocate matrix1\n");
        exit(1);
    }
	
// aloca o eixo final, ataves de um apontador de arrays
    for(x = 0; x < N; x++) {
        matrix.data[x] = (char **) malloc(N * sizeof(char *));
        if(matrix.data[x] == NULL) {
            printf("Failed to allocate matrix3\n");
            exit(1);
        }
        for (y = 0; y < N; y++){
            matrix.data[x][y] = (char*) calloc(N, sizeof(char));
            if(matrix.data[x][y] == NULL) {
                printf("Failed to allocate matrix6\n");
                exit(1);
            }
            for (z = 0; z < N; z++)
                if(r4_uni() < density)
                    {
						// preenchimento initial do grid_even dependendo da seed
                        matrix.data[x][y][z] = (int)(r4_uni() * N_SPECIES) + 1; // preenchimento initial do grid_even dependendo da seed
                    }
        }     
    }
}

void flattenMatrix(Matrix3D *matrix, char *flatMatrix, int Z_DIM) {
    for (int i = 0; i < X_DIM; i++) {
        for (int j = 0; j < Y_DIM; j++) {
            for (int k = 0; k < Z_DIM; k++) {
                flatMatrix[i * Y_DIM * Z_DIM + j * Z_DIM + k] = matrix->data[i][j][k];
            }
        }
    }
}

void reconstructMatrix(char *flatMatrix, char received_matrix[X_DIM][Y_DIM]) {
    for (int i = 0; i < X_DIM; i++) {
        for (int j = 0; j < Y_DIM; j++) {

            received_matrix[i][j] = flatMatrix[i * Y_DIM + j];
            
        }
    }
}

void printMatrix(char received_matrix[X_DIM][Y_DIM], int rank) {
    for (int i = 0; i < X_DIM; i++) {
        for (int j = 0; j < Y_DIM; j++) {
           
            printf("Process %d: matrix[%d][%d] = %d\n", rank, i, j, received_matrix[i][j]);
        }
    }
    printf("\n");
}

int main(int argc, char *argv[]) {

    int number_of_gens,  number_of_cells;
    float density;

    number_of_gens = atoi (argv[1]);
    number_of_cells = atoi (argv[2]);
    density = atof (argv[3]);
    seed = atoi (argv[4]);

    printf("gens: %d    cells:%d    density:%f  seed:%d\n",number_of_gens, number_of_cells, density, seed );

    int Z_DIM = number_of_cells;
    
    int rank, size;
    int sendcounts[X_DIM], displs[X_DIM];
    //char flatMatrix[X_DIM * Y_DIM * Z_DIM];
    matrix.data = gen_initial_grid(number_of_cells, density, seed);
    
    char * flatMatrix = (char *) malloc((X_DIM*Y_DIM*Z_DIM) * sizeof(char));
        if(flatMatrix == NULL) {
            printf("Failed to allocate matrix3\n");
            exit(1);
        }

    for(int aux_x = 0; aux_x < (X_DIM*Y_DIM*Z_DIM); aux_x++) {
        
        flatMatrix[aux_x]=0;
    }
    

    char recvBuffer[X_DIM * Y_DIM]={0}; // Separate receive buffer for each process
    char received_matrix[X_DIM][Y_DIM];
    

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Datatype matrix3d_type;
    MPI_Type_contiguous(X_DIM * Y_DIM * Z_DIM, MPI_INT, &matrix3d_type);
    MPI_Type_commit(&matrix3d_type);

    

    for (int i = 0; i < X_DIM; i++) {
        for (int j = 0; j < Y_DIM; j++) {

                received_matrix[i][j] = 0;
   
        }
    }

    // Prepare data to scatter
    if (rank == 0) {
        // Initialize the matrix
        for (int i = 0; i < X_DIM; i++) {
            for (int j = 0; j < Y_DIM; j++) {
                for (int k = 0; k < Z_DIM; k++) {
                    matrix.data[i][j][k] = (i * Y_DIM * Z_DIM) + (j * Z_DIM) + k; // Ensure values are within manageable range
                }
            }
        }

        // Print the initial matrix
        printf("Initial matrix (rank 0):\n");
            for(int aux_x4=0; aux_x4<X_DIM; aux_x4++)
            {
                for(int aux_y4=0; aux_y4<Y_DIM; aux_y4++)
                {
                    for(int aux_z4=0; aux_z4<Z_DIM; aux_z4++)
                    {
                        printf("%d ", matrix.data[aux_x4][aux_y4][aux_z4]);
                    }
                    printf("\n");
                }
                printf("\n");
            }

        // Flatten the matrix
        flattenMatrix(&matrix, flatMatrix, Z_DIM);

        // Prepare sendcounts and displs for scatterv
        for (int i = 0; i < X_DIM; i++) {
            sendcounts[i] = Y_DIM * Z_DIM;
            displs[i] = i * Y_DIM * Z_DIM;
        }
        
        // Print sendcounts and displs for debugging
        printf("Send counts (rank 0): ");
        for (int i = 0; i < X_DIM; i++) {
            printf("%d ", sendcounts[i]);
        }
        printf("\n");
        printf("Displacements (rank 0): ");
        for (int i = 0; i < X_DIM; i++) {
            printf("%d ", displs[i]);
        }
        printf("\n");
    }

    // Ensure all processes reach this point before proceeding
    MPI_Barrier(MPI_COMM_WORLD);

    // Scatter the matrix to all processes
    MPI_Scatterv(flatMatrix, sendcounts, displs, MPI_INT, recvBuffer, Y_DIM * Z_DIM, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    // Print received data for debugging
    printf("Received data (rank %d): ", rank);
    for (int i = 0; i < Y_DIM * Z_DIM; i++) {
        printf("%d ", recvBuffer[i]);
    }
    printf("\n");

    // Ensure all processes receive their data before proceeding
    MPI_Barrier(MPI_COMM_WORLD);

    // Reconstruct the matrix
    reconstructMatrix(recvBuffer, received_matrix);

    MPI_Barrier(MPI_COMM_WORLD);
    // Print the received matrix for each process
    printf("Received matrix (rank %d):\n", rank);
    printMatrix(received_matrix, rank);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Type_free(&matrix3d_type);
    MPI_Finalize();
    return 0;
}
