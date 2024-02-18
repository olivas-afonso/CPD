#include <stdio.h>
#include <stdlib.h>

#define N_SPECIES 9

char *** grid;
unsigned int seed;
long long cont_species[9]={0,0,0,0,0,0,0,0,0};
//long long cont1=0, cont2=0, cont3=0, cont4=0, cont5=0, cont6=0, cont7=0, cont8=0, cont9=0;
//long long max_cont1, max_cont2, max_cont3, max_cont4, max_cont5, max_cont6, max_cont7, max_cont8, max_cont9;
long long max_count[9]={0,0,0,0,0,0,0,0,0};
int max_gen[9];
int current_gen =0;
//int max_gen1, max_gen2, max_gen3, max_gen4, max_gen5, max_gen6, max_gen7, max_gen8, max_gen9; 

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

void search_species_line(long long N, char * line) // search linha
{
    long long aux;

    for(aux =0; aux<N; aux++)
    {
        switch (line[aux])
        {
        case 1:
            cont_species[0]++;
            break;
        
        case 2:
            cont_species[1]++;
            break;

        case 3:
            cont_species[2]++;
            break;

        case 4:
            cont_species[3]++;
            break;

        case 5:
            cont_species[4]++;
            break;

        case 6:
            cont_species[5]++;
            break;

        case 7:
            cont_species[6]++;
            break;

        case 8:
            cont_species[7]++;
            break;
        
        case 9:
            cont_species[8]++;
            break;
        
        default:
            break;
        }
    }
    
}

void life_rules(long long N, char *** grid)
{
    long long aux;
    for(aux=0; aux< N; aux ++)
    {
        
    }
}

int main(int argc, char *argv[]) {
   // printf() displays the string inside quotation
    long long contx, conty; 
    int auxi;
   // primeiro numero gen, segundo nº celulas por lado do cubo, terceiro densidade initial, quarto seed
   printf("Hello, World!\n");
   //printf("%s\n", argv[1]);

   grid = gen_initial_grid(4, .4, 100);
   printf("Value :%d\n", grid[0][0][0]);

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

        printf("%d %llu %d\n", auxi+1, max_count[auxi], max_gen[auxi]);
    }

    
    //printf("Count :%llu\n", cont_species[2]);

   return 0;
}
