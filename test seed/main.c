#include <stdio.h>
#include <stdlib.h>

#include "soma.h" // Inclui o arquivo de header para usar a função soma


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

char ***gen_initial_grid(long long N, float density, int input_seed)
{
    int x, y, z;
    
    grid = (char ***) malloc(N * sizeof(char **));
    if(grid == NULL) {
        printf("Failed to allocate matrix\n");
        exit(1);
    }
    for(x = 0; x < N; x++) {
        grid[x] = (char **) malloc(N * sizeof(char *));
        if(grid[x] == NULL) {
            printf("Failed to allocate matrix\n");
            exit(1);
        }
        grid[x][0] = (char *) calloc(N * N, sizeof(char)); // N by N cube 
        if(grid[x][0] == NULL) {
            printf("Failed to allocate matrix\n");
            exit(1);
        }
        for (y = 1; y < N; y++)
            grid[x][y] = grid[x][0] + y * N;
    }

    init_r4uni(input_seed);
    for (x = 0; x < N; x++)
        for (y = 0; y < N; y++)
            for (z = 0; z < N; z++)
                if(r4_uni() < density)
                    grid[x][y][z] = (int)(r4_uni() * N_SPECIES) + 1; // preenchimento initial do grid dependendo da seed

    return grid;
}


int main(int argc, char *argv[]) {
   // printf() displays the string inside quotation
    long long contx, conty; 
    int auxi;
	//int ger;
   // primeiro numero gen, segundo nº celulas por lado do cubo, terceiro densidade initial, quarto seed
   printf("jogo da vida!\n");
   fflush(stdout);
   
   grid = gen_initial_grid(128, .5, 1000);
	
	
//	ver_grid(4);
  
  
	//ver_grid(4);
	for(int n_gera = 0; n_gera < 201; n_gera++){
		printf("geracao:%d\n", n_gera);
		geracao();
	}
   //printf("ACABOU A GERACAO\n\n");
   
  // ver_grid(4);
   
   //ver_vizinhos();
	/*
    for(contx=0; contx < 4; contx++){
        for(conty=0; conty < 4; conty++){
            search_species_line(4, grid[contx][conty]);
        }
    }

    for(auxi=0; auxi < 9; auxi++)
    {
        if(cont_species[auxi] > max_count[auxi])
        {
           max_count[auxi] = cont_species[auxi];
           max_gen[auxi]=current_gen;
        }
	
   //     printf("%d %llu %d\n", auxi+1, max_count[auxi], max_gen[auxi]);
    }
	*/
    
    //printf("Count :%llu\n", cont_species[2]);

   return 0;
}
