/*for(int merdoca = 0; merdoca < 64; merdoca++){
					aux++;
					//printf("celula: %d",aux);
					// preenchimento initial do grid_even dependendo da seed
					prov = (int)(r4_uni() * N_SPECIES) + 1; // preenchimento initial do grid_even dependendo da seed
					//count_species[grid_even[x][y][z]]++;
					
					if(merdoca < 2){
						grid_even[x][y][z] = prov;
						x++;	
					}
					y++;
					
					if(merdoca>3 && merdoca < 6){
						grid_even[x][y][z] = prov;
					}
					z++;
					x = 0;
					y = 0;
					
					if(merdoca>15 && merdoca < 18){
						grid_even[x][y][z] = prov;
					}
					x++;
					if(merdoca>19 && merdoca < 23){
						grid_even[x][y][z] = prov;
					}
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
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
						printf("numero: %d\n",grid_even[x][y][z]);
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


int main(int argc, char *argv[]) {
    double exec_time;
    int auxi;
    int gen_number=0;

    int aux_x4, aux_y4, aux_z4; 

    int number_of_gens,  number_of_cells;
    float density;
    
    // retira os argumentos do terminal e colaca em variaveis
    number_of_gens = atoi (argv[1]);
    number_of_cells = atoi (argv[2]);
    density = atof (argv[3]);
    seed = atoi (argv[4]);

    //cria a grid aleatoria atraves dos inputs (funcao fornecida)
    grid_even = gen_initial_grid(number_of_cells, density, seed);
  
  return 0;
}