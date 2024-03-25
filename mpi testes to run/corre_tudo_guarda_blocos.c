#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <omp.h>
#include <mpi.h>

#define N_SPECIES 9


char *** grid_even;
char *** grid_odd;
unsigned int seed;
long count_species[10]={0,0,0,0,0,0,0,0,0,0};

long max_count[10]={0,0,0,0,0,0,0,0,0};
int max_gen[10];


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


/************************************************************************************************
* Nome: gen_initial_grid
* funcao: gera as duas matrizes, a inicial e a auxiliar; ao longo das geracoes as celulas iram
*passar de uma para a outra.
************************************************************************************************/
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
                        count_species[grid_even[x][y][z]]++;
						
						//preenchimento_da_grid(grid_even[x][y][z]);
						
					
					
					}
        }     
    }

	// conta as especies da geracao 0
    for(x=1; x < 10; x++)
    {
        if(count_species[x] > max_count[x])
        {
            max_count[x] = count_species[x];
            max_gen[x]=0;
        }
    } 
                    

    return grid_even;
}

/************************************************************************************************
* Nome: gen_initial_bloco
* funcao: gera as duas matrizes, a inicial e a auxiliar; ao longo das geracoes as celulas iram
*passar de uma para a outra.
************************************************************************************************/
char ***gen_initial_bloco(int max, int min ,int n,int N, float density, int input_seed)
{
    int x, y, z;
	int prov = 0;
    int aux = 0;
	
    init_r4uni(input_seed);
    
//alocacao da memeoria dinamica, alocando primeiro um apontador triplo o que corresponde a uma dimensao do cubo  ~
    grid_even = (char ***) malloc(max* sizeof(char **));
    if(grid_even == NULL) {
        printf("Failed to allocate matrix1\n");
        exit(1);
    }
	
//aloca a dimensao x atraves de um apontador de um apontador    
    /*grid_odd = (char ***) malloc(N * sizeof(char **));
    if(grid_odd == NULL) {
        printf("Failed to allocate matrix2\n");
        exit(1);
    }*/
// aloca o eixo final, ataves de um apontador de arrays
    for(x = 0; x < n; x++) {
        grid_even[x] = (char **) malloc(n * sizeof(char *));
        if(grid_even[x] == NULL) {
            printf("Failed to allocate matrix3\n");
            exit(1);
        }

        /*grid_odd[x] = (char **) malloc(N * sizeof(char *));
        if(grid_odd[x] == NULL) {
            printf("Failed to allocate matrix5\n");
            exit(1);
        }*/

        for (y = 0; y < n; y++){
            grid_even[x][y] = (char*) calloc(n, sizeof(char));
            if(grid_even[x][y] == NULL) {
                printf("Failed to allocate matrix6\n");
                exit(1);
            }
          /*  grid_odd[x][y] = (char*) calloc(N, sizeof(char));
            if(grid_odd[x][y] == NULL) {
                printf("Failed to allocate matrix6\n");
                exit(1);
            }*/
            for (z = 0; z < n; z++)
				grid_even[x][y][z] = 0;
        }     
    }

	x = 0;
	y = 0;
	z = 0;	

/*
	for(aux = 0; aux <64; aux++){
		if(r4_uni() < density){
			prov = (int)(r4_uni() * N_SPECIES) + 1; // preenchimento initial do grid_even dependendo da seed
			//count_species[grid_even[x][y][z]]++;
			printf("numero: %d\n",prov);	
		}	
		if(aux < 2){
				grid_even[x][y][z] = prov;
				x++;	
			}
			y++;
			
			if(aux>3 && aux < 6){
				grid_even[x][y][z] = prov;
			}
			z++;
			x = 0;
			y = 0;
			
			if(aux>15 && aux < 18){
				grid_even[x][y][z] = prov;
			}
			x++;
			if(aux>19 && aux < 23){
				grid_even[x][y][z] = prov;
			}
	}
*/
	for(x = 0; x < N; x++){
		for(y=0; y < N; y++){
			for(z=0; z<N; z++){
				if(r4_uni() < density){
					prov = (int)(r4_uni() * N_SPECIES) + 1; // preenchimento initial do grid_even dependendo da seed
					//count_species[grid_even[x][y][z]]++;
					if (x>= min && y>= min && z>= min) {
						if (x<max && y<max && z<max) {
							grid_even[x][y][z] = prov;
							printf("numero: %d\n",prov);
						}	
					}
				}		
			}
		}
		
	}



	// conta as especies da geracao 0
    for(x=1; x < 10; x++)
    {
        if(count_species[x] > max_count[x])
        {
            max_count[x] = count_species[x];
            max_gen[x]=0;
        }
    } 
                    

    return grid_even;
}

/************************************************************************************************
* Nome:death_rule
* funcao: verifica os vizinhos no caso da celula ter estado morta na geracao anterior
*
************************************************************************************************/
int death_rule(int N, char *** grid, long aux_x, long aux_y, long aux_z)
{

    long search_x, search_y, search_z;
    long aux_search_y, aux_search_z;
    long cont_species_death[9]={0,0,0,0,0,0,0,0,0};
    int cont_rule=0;
    int max=0, max_pos=0, i;
    int x,y,z;
    
    for(search_x= (aux_x-1+N)%N, x=0; x < 3; x++, search_x++) 
    {
        for(search_y=(aux_y-1+N)%N, y=0; y < 3; y++, search_y++)
        {
            for(search_z=(aux_z-1+N)%N, z=0; z< 3;z++, search_z++)
            {
                if (grid [search_x % N][search_y % N][search_z % N] != 0){       
                    ++cont_rule;                
                    cont_species_death[grid[search_x % N][search_y% N][search_z% N]-1]++;
                }
                
                if (cont_rule >10){
                    return 0;
                }               
            }
        }
    } 
    
     if ( cont_rule >= 7 && cont_rule <= 10 )
    {
        max=cont_species_death[0];
        
        max_pos=0;
        for(i=1; i < 9;i++ )
        {
            if(cont_species_death[i]>max)
            {
                max = cont_species_death[i];
                max_pos=i;         
            }      
        }

        return max_pos+1;
    }
    else return 0; 
}

/************************************************************************************************
* Nome:life_rule
* funcao: verifica os vizinhos no caso da celula ter estado viva na geracao anterior
*
************************************************************************************************/
int life_rule (int N, char *** grid, long aux_x, long aux_y, long aux_z){
    long search_x, search_y, search_z;
    long aux_search_y, aux_search_z;
    int cont_rule=-1;
    int x,y,z;
	
	// corre vizinhos e em caso de extremo verifica o extremo oposto
    for(search_x= (aux_x-1+N)%N, x=0; x < 3; x++, search_x++) 
    {
        for(search_y=(aux_y-1+N)%N, y=0; y < 3; y++, search_y++)
        {
            for(search_z=(aux_z-1+N)%N, z=0; z< 3;z++, search_z++)
            {
				//verifica se o vizinho está vivo
                if (grid [search_x % N][search_y % N][search_z % N] != 0){       
                    ++cont_rule;
                }
				// (OTIMIZACAO) se já tem vizinhos suficientes para permanecer morta sai da funcao
                if (cont_rule >13){
                    return 0;
                }               
            }
        }
    } 
	 
	// se nao tiver celulas suficientes permanece morto
    if (cont_rule<=4){
        return 0;
    }else{
		// revive a celula	
        return grid[aux_x][aux_y][aux_z];
    }
}

/************************************************************************************************
* Nome:rules
* funcao: corre as funcoes relativas às regras dependendo do valor que a celula tem
*
************************************************************************************************/
void rules(int N, char ***grid_new, char ***grid_old)
{
    long aux_x, aux_y, aux_z;

	// corre todas as celulas da grid 
    for(aux_x=0; aux_x< N; aux_x ++)
    {
        for(aux_y=0; aux_y<N; aux_y++)
        {
            for(aux_z=0; aux_z<N; aux_z++)
            {
                if(grid_old[aux_x][aux_y][aux_z]==0) // morto 
                { 
                    grid_new[aux_x][aux_y][aux_z]= death_rule(N, grid_old, aux_x, aux_y, aux_z);
                }
                else
                {  
                    grid_new[aux_x][aux_y][aux_z]= life_rule(N, grid_old, aux_x, aux_y, aux_z);     
                }
				// se a celula esta viva nesta geracao, aumentamos o numero no array contador 
                count_species[grid_new[aux_x][aux_y][aux_z]]++;
            }
        }
    }
}

/************************************************************************************************
* Nome:liberar matriz
* funcao: liberta a memoria da grid auxiliar 
*
************************************************************************************************/
void freeMatrix(int N) {
    int i, j;

	
    // Libera a memória da terceira dimensão
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            free(grid_even[i][j]);
            free(grid_odd[i][j]);
        }
        free(grid_even[i]);
        free(grid_odd[i]);
    }
	// Libera a memória da primeira dimensão
    free(grid_even);
    free(grid_odd);
}






double raiz_cubica_funcao(double num) {
    double x0 = num / 3; // Estimativa inicial
    double x1 = (2 * x0 + num / (x0 * x0)) / 3; // Melhor estimativa

    // Continuar refinando até convergência
    while (fabs(x1 - x0) >= 0.00001) {
        x0 = x1;
        x1 = (2 * x0 + num / (x0 * x0)) / 3;
    }

    return x1;
}

/************************************************************************************************
* Nome:main
* funcao: corre as funcoes pela ordem de execução correta
*
************************************************************************************************/
int main(int argc, char *argv[]) {
    
	/*************************************************************************
	* Var do serial
	***************************************************************************/
	double exec_time;
    int auxi;		// contador do numero de especies por geração
    int gen_number=0;

    int aux_x4, aux_y4, aux_z4; 

    int number_of_gens,  number_of_cells;
    float density;


	long aux_x, aux_y, aux_z;


	/*************************************************************************
	* Var do MPI
	***************************************************************************/	
	int rank, size;
    int dims[3] = {0, 0, 0};
    int periods[3] = {1, 1, 1};
    int coords[3];
    MPI_Comm cart_comm;
    MPI_Status status;
	
	
	// retira os argumentos do terminal e colaca em variaveis
    number_of_gens = atoi (argv[1]);
    number_of_cells = atoi (argv[2]);
    density = atof (argv[3]);
    seed = atoi (argv[4]);
	
	
	double N_double;
	double n_double;
	double x_double;
    
	
	int NX = number_of_cells;
	int NY = NX;
	int NZ = NX;
	
	printf("valor da aresta = %d\n", NZ);
	
	
    // Create a 3D array to hold the layer of the grid for each process
	//int layer[NX][NY][NZ];
	
	/*************************************************************************
	* codigo do MPI
	***************************************************************************/	

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
    

	
	/*************************************************************************
	* codigo do serial 
	***************************************************************************/
	
if(rank == 1){
	x_double = (double)size;
	N_double = (double)NX;
	
	printf("Você digitou: %.2lf e %.2lf\n", N_double, x_double);
	

	printf("x = %2.f \n", x_double);
	double raiz_cubica_x = raiz_cubica_funcao(x_double);

	printf("raiz = %2.f \n", raiz_cubica_x);

    // Calcular n
	//double = (double)N;
	printf("N = %2.f \n", N_double);
	
    n_double = N_double / raiz_cubica_x;

    // Exibir o resultado
    printf("O tamanho dos cubos pequenos (n) é: %.2f %.2f\n ", n_double, N_double);
	
	double comprimento_cubo_pequeno = N_double / raiz_cubica_x;
	
	printf("valor de n=%f\n", comprimento_cubo_pequeno);

}

	
	
	//cria a grid aleatoria atraves dos inputs (funcao fornecida)
for (int i = 0; i < size; i++) {
	if (rank == 1) {
		char *** grid;
		int min = 0;
		int max = 0;
		int inicio = 0;
		int n = 0;
		
		n = (int)n_double;
		min = ((int)n_double)*rank;
		max = min + (((int)n_double));

/*		
		inicio = ((int)n_double)*rank; 
		
		*x = inicio % 4;
		*y = (inicio / 4) % 4;
		*z = inicio / 16;
		
		min = x + y * 4 + z * 16;
		max = x + y * 4 + z * 16; 
*/
		printf("valores: %d %d %d\n", min ,max, n);	
		
		grid = gen_initial_bloco(max, min,n,number_of_cells, density, seed);
		
		
		
		printf("Rank %d: Layer %d\n", rank, rank);	
			for(int merda_z= 0; merda_z < max; merda_z++){
				for(int merda_y = 0; merda_y < max; merda_y++){
					for(int merda_x = 0; merda_x < max; merda_x++){
						printf("%d ", grid[merda_z][merda_y][merda_x]);
					}			
						printf("\n");			
				}
				printf("\n\n");
			}
		
		
		
		
		//exec_time = -omp_get_wtime();

		/*// corre todas as celulas da grid 
		for(aux_x=0; aux_x< NX; aux_x ++)
		{
			for(aux_y=0; aux_y<NY; aux_y++)
			{
				for(aux_z=0; aux_z<NZ; aux_z++)
				{
					 layer[aux_x][aux_y][aux_z] = grid[aux_x][aux_y][aux_z]; 
				}
			}
		
		}
			printf("Rank %d: Layer %d\n", rank, rank);	
			for(int merda_z= 0; merda_z < NX; merda_z++){
				for(int merda_y = 0; merda_y < NX; merda_y++){
					for(int merda_x = 0; merda_x < NX; merda_x++){
						printf("%d ", grid[merda_z][merda_y][merda_x]);
					}			
						printf("\n");			
				}
				printf("\n\n");
			}
		
		
		freeMatrix(16);
		MPI_Barrier(MPI_COMM_WORLD);
	*/
	}  
}	
	
	
   // exec_time += omp_get_wtime();
    //fprintf(stderr, "%.1fs\n", exec_time);
    
    // Synchronize the output
    MPI_Barrier(MPI_COMM_WORLD);
/*
    // Only let process 0 print the initial grid
    if (rank == 0) {
        printf("Initial Grid:\n");
        for (int k = 0; k < NZ; k++) {
			printf("Rank %d: Layer %d\n", rank, rank);
				for (int i = 0; i < NY; i++) {
					for (int j = 0; j < NZ; j++) {
						printf("%2d ", layer[k][i][j]);
					}
					printf("\n");
				}
				printf("\n");
				
		}
    }

    // Synchronize the output
    MPI_Barrier(MPI_COMM_WORLD);
/*
    // Print out the layer of the grid for each process one at a time
    for (int i = 0; i < size; i++) {
        if (rank == i) {
            printf("Rank %d is printing...\n", rank);
            for (int k = start_layer; k < end_layer; k++) {
                MPI_Barrier(MPI_COMM_WORLD); // Synchronize before receiving
                MPI_Bcast(&layer, NX * NY * NZ, MPI_INT, i, MPI_COMM_WORLD); // Broadcast the layer
                
				printf("Rank %d: Layer %d\n", rank, rank);
				for (int i = 0; i < NZ; i++) {
					for (int j = 0; j < NY; j++) {
						printf("%2d ", layer[i][j][rank]);
					}
					printf("\n");
				}
				printf("\n");
				
				
				if (k != end_layer - 1) {
                    printf("\n"); // Print new line unless it's the last layer
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
*/
    MPI_Comm_free(&cart_comm);
    MPI_Finalize();
	

    return 0;
}


