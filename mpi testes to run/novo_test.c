#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>


#define X_DIM 3
#define Y_DIM 3
#define Z_DIM 3

#define N_SPECIES 9


char *** grid_even;
char *** grid_odd;
unsigned int seed;

typedef struct {
    char data[X_DIM][Y_DIM][Z_DIM];
} Matrix3D;

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

char ***gen_initial_grid(int N, float density, int input_seed)
{
    int x, y, z;
    
    init_r4uni(input_seed);
    
//alocacao da memeoria dinamica, alocando primeiro um apontador triplo o que corresponde a uma dimensao do cubo  ~
    grid_even = (char ***) malloc(N * sizeof(char **));
    if(grid_even == NULL) {
        printf("Failed to allocate matrix1\n");
        exit(1);
    }
	
//aloca a dimensao x atraves de um apontador de um apontador    
    grid_odd = (char ***) malloc(N * sizeof(char **));
    if(grid_odd == NULL) {
        printf("Failed to allocate matrix2\n");
        exit(1);
    }
// aloca o eixo final, ataves de um apontador de arrays
    for(x = 0; x < N; x++) {
        grid_even[x] = (char **) malloc(N * sizeof(char *));
        if(grid_even[x] == NULL) {
            printf("Failed to allocate matrix3\n");
            exit(1);
        }

        grid_odd[x] = (char **) malloc(N * sizeof(char *));
        if(grid_odd[x] == NULL) {
            printf("Failed to allocate matrix5\n");
            exit(1);
        }

        for (y = 0; y < N; y++){
            grid_even[x][y] = (char*) calloc(N, sizeof(char));
            if(grid_even[x][y] == NULL) {
                printf("Failed to allocate matrix6\n");
                exit(1);
            }
            grid_odd[x][y] = (char*) calloc(N, sizeof(char));
            if(grid_odd[x][y] == NULL) {
                printf("Failed to allocate matrix6\n");
                exit(1);
            }
            for (z = 0; z < N; z++)
                if(r4_uni() < density)
                    {
						// preenchimento initial do grid_even dependendo da seed
                        grid_even[x][y][z] = (int)(r4_uni() * N_SPECIES) + 1; // preenchimento initial do grid_even dependendo da seed
                    }
        }     
    }
    /*
	// conta as especies da geracao 0
    for(x=1; x < 10; x++)
    {
        if(count_species[x] > max_count[x])
        {
            max_count[x] = count_species[x];
            max_gen[x]=0;
        }
    } 
    */

   
                    

    return grid_even;
}

int main(int argc, char *argv[])
{
    int rank, size;
    int sendcounts, displs;
    MPI_Comm comm;
    Matrix3D matrix, received_matrix;
    
    int dim[2], period[2], reorder;
    int coord[2], id;
    int number_of_gens,  number_of_cells;
    float density;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    number_of_gens = atoi (argv[1]);
    number_of_cells = atoi (argv[2]);
    density = atof (argv[3]);
    seed = atoi (argv[4]);

    grid_even = gen_initial_grid(number_of_cells, density, seed);
    /*
    char send_stuff[NUM_LINES][NUM_COL][number_of_cells];

    int i, j, k;
    for (i = 0; i < NUM_LINES; i++) {
        for (j = 0; j < NUM_COL; j++) {
            for (k = 0; k < number_of_cells; k++) {
                send_stuff[i][j][k] = 0; // Fill each element with zero
            }
        }
    }
    */

    /*
    MPI_Type_vector(NUM_LINES, NUM_COL, 4, MPI_CHAR, &blocktype2);
    MPI_Type_create_resized( blocktype2, 0, sizeof(char), &blocktype);
    MPI_Type_commit(&blocktype);
    */

    MPI_Datatype matrix3d_type;
    MPI_Type_contiguous(X_DIM * Y_DIM * Z_DIM, MPI_INT, &matrix3d_type);
    MPI_Type_commit(&matrix3d_type);

    if (rank == 0) {
        for (int i = 0; i < X_DIM; i++) {
            for (int j = 0; j < Y_DIM; j++) {
                for (int k = 0; k < Z_DIM; k++) {
                    matrix.data[i][j][k] = rank + i + j + k;
                }
            }
        }
        sendcounts = X_DIM * Y_DIM * Z_DIM;
        displs = 0;
    }

    MPI_Scatterv(grid_even, &sendcounts, &displs, matrix3d_type, &received_matrix, 1, matrix3d_type, 0, MPI_COMM_WORLD);

    // Each process prints the received data
    for (int x = 0; x < X_DIM; x++) {
        for (int y = 0; y < Y_DIM; y++) {
            for (int z = 0; z < Z_DIM; z++) {
                printf("Process %d: received matrix[%d][%d][%d] = %d\n", rank, x, y, z, received_matrix.data[x][y][z]);
            }
        }
    }

    MPI_Type_free(&matrix3d_type);
    MPI_Finalize();


}