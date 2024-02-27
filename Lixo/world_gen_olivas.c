#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N_SPECIES 9

char *** grid_even;
char *** grid_odd;
unsigned int seed;
long count_species[9]={0,0,0,0,0,0,0,0,0};

//long long cont1=0, cont2=0, cont3=0, cont4=0, cont5=0, cont6=0, cont7=0, cont8=0, cont9=0;
//long long max_cont1, max_cont2, max_cont3, max_cont4, max_cont5, max_cont6, max_cont7, max_cont8, max_cont9;
long max_count[9]={0,0,0,0,0,0,0,0,0};
long max_gen[9];

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
    
    grid_even = (char ***) malloc(N * sizeof(char **));
    if(grid_even == NULL) {
        printf("Failed to allocate matrix\n");
        exit(1);
    }
    for(x = 0; x < N; x++) {
        grid_even[x] = (char **) malloc(N * sizeof(char *));
        if(grid_even[x] == NULL) {
            printf("Failed to allocate matrix\n");
            exit(1);
        }
        grid_even[x][0] = (char *) calloc(N * N, sizeof(char)); // N by N cube 
        if(grid_even[x][0] == NULL) {
            printf("Failed to allocate matrix\n");
            exit(1);
        }

        for (z=0)
        for (y = 1; y < N; y++)
            grid_even[x][y] = grid_even[x][0] + y * N;
    }

    init_r4uni(input_seed);
    for (x = 0; x < N; x++)
        for (y = 0; y < N; y++)
            for (z = 0; z < N; z++)
                if(r4_uni() < density)
                    grid_even[x][y][z] = (int)(r4_uni() * N_SPECIES) + 1; // preenchimento initial do grid_even dependendo da seed

    
    
    
    grid_odd = (char ***) malloc(N * sizeof(char **));
    if(grid_odd == NULL) {
        printf("Failed to allocate matrix\n");
        exit(1);
    }
    for(x = 0; x < N; x++) {
        grid_odd[x] = (char **) malloc(N * sizeof(char *));
        if(grid_odd[x] == NULL) {
            printf("Failed to allocate matrix\n");
            exit(1);
        }
        grid_odd[x][0] = (char *) calloc(N * N, sizeof(char)); // N by N cube 
        if(grid_odd[x][0] == NULL) {
            printf("Failed to allocate matrix\n");
            exit(1);
        }
        for (y = 1; y < N; y++)
            grid_odd[x][y] = grid_odd[x][0] + y * N;
    }
    return grid_even;

}

void liberarMatriz(int N) {
    int i, j;

    // Libera a memória da terceira dimensão
    for (i = 0; i < N; i++) {
            free(grid_even[i][0]);
            free(grid_odd[i][0]);
    }

    // Libera a memória da segunda dimensão
    for (i = 0; i < N; i++) {
        free(grid_even[i]);
        free(grid_odd[i]);
    }

    // Libera a memória da primeira dimensão
    free(grid_even);
    free(grid_odd);
}


void count_gen0(long long N, char *** grid) // search linha
{
    int aux_x, aux_y, aux_z, auxi;

    for(aux_x=0; aux_x<N; aux_x++)
    {
        for(aux_y=0; aux_y<N; aux_y++)
        {
            for(aux_z=0; aux_z<N; aux_z++)
            {   
                if(grid[aux_x][aux_y][aux_z] !=0)
                count_species[grid[aux_x][aux_y][aux_z]-1]++;
            }
        }
    }
    for(auxi=0; auxi < 9; auxi++)
        {
            if(count_species[auxi] > max_count[auxi])
            {
                max_count[auxi] = count_species[auxi];
                max_gen[auxi]=0;
            }
        }
    
}

int death_rule(long long N, char *** grid, long aux_x, long aux_y, long aux_z)
{

    long search_x, search_y, search_z;
    long aux_search_y, aux_search_z;
    long cont_species_death[9]={0,0,0,0,0,0,0,0,0};
    int cont_rule=0;
    int max=0, max_pos=0, i;
    int x,y,z;
    
    if (aux_x == 0)
    {
        search_x = N-1;
    }
    else{
        search_x = aux_x-1;  
    }

    if (aux_y == 0)
    {    
        search_y = N-1;
        aux_search_y = search_y;
    }
    else{
        search_y = aux_y-1;
        aux_search_y = search_y; 
    }

    if (aux_z == 0)
    {
        search_z = N-1;
        aux_search_z = search_z;
    }
    else{
        
        search_z = aux_z-1;
        aux_search_z = search_z; 
    }

    for(x=0; x < 3; x++, search_x++) 
    {
        for(search_y=aux_search_y, y=0; y < 3; y++, search_y++)
        {
            for(search_z=aux_search_z, z=0; z< 3;z++, search_z++)
            {
                if ((search_x % N) != aux_x || (search_y % N) != aux_y || (search_z % N) != aux_z){
                    if (grid [search_x % N][search_y % N][search_z % N] != 0){       
                        ++cont_rule;
                        cont_species_death[grid[search_x % N][search_y% N][search_z% N]-1]++;
                         
                    }
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


int life_rule (long long N, char *** grid, long aux_x, long aux_y, long aux_z){
    long search_x, search_y, search_z;
    long aux_search_y, aux_search_z;
    int cont_rule=0;
    int x,y,z;
    
    if (aux_x == 0)
    {
        search_x = N-1;
    }
    else{
        search_x = aux_x-1;  
    }

    if (aux_y == 0)
    {
        search_y = N-1;
        aux_search_y = search_y;
    }
    else{
        search_y = aux_y-1;
        aux_search_y = search_y; 
    }

    if (aux_z == 0)
    {
        search_z = N-1;
        aux_search_z = search_z;
    }
    else{
        
        search_z = aux_z-1;
        aux_search_z = search_z; 
    }

   
    for(x=0; x < 3; x++, search_x++) 
    {
        for(search_y=aux_search_y, y=0; y < 3; y++, search_y++)
        {
            for(search_z=aux_search_z, z=0; z< 3;z++, search_z++)
            {
                if ((search_x % N) != aux_x || (search_y % N) != aux_y || (search_z % N) != aux_z){
                    if (grid [search_x % N][search_y % N][search_z % N] != 0){       
                        ++cont_rule;
                    }
                }
                if (cont_rule >13){
                    return 0;
                }               
            }
        }
    } 

    if (cont_rule<=4){
        return 0;
    }else{  
        return grid[aux_x][aux_y][aux_z];
    }
}

void rules(long long N, char ***grid_new, char ***grid_old)
{
    long aux_x, aux_y, aux_z;
    long long cont_rule=0;

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

                if (grid_new[aux_x][aux_y][aux_z] != 0){
                    //printf("x:%d y:%d z:%d  value:%d\n", aux_x, aux_y, aux_z, grid_new[aux_x][aux_y][aux_z]);
                    count_species[grid_new[aux_x][aux_y][aux_z]-1]++;
                }
            }
        }
    }
}


int main(int argc, char *argv[]) {
   // printf() displays the string inside quotation
    int auxi;
    int gen_number=0;

    //int aux_x4, aux_y4, aux_z4; 

    long number_of_gens, number_of_cells;
    float density;
    
    number_of_gens = atoi (argv[1]);
    number_of_cells = atoi (argv[2]);
    density = atof (argv[3]);
    seed = atoi (argv[4]);
    printf("%d %d %f %d \n \n", number_of_gens ,number_of_cells, density, seed);
    grid_even = gen_initial_grid(number_of_cells, density, seed);
    count_gen0(number_of_cells, grid_even);
    //printf("1: %d current_gen: %d\n", count_species[0], gen_number);

    for(gen_number=1; gen_number<=number_of_gens; gen_number++)
    {
        for(auxi=0; auxi < 9; auxi++)
        {
            count_species[auxi]=0;
        }
        
        if(gen_number %2 == 1) // write to odd
        {
        
            rules(number_of_cells, grid_odd, grid_even);
            /*
            printf("GENERATION %d\n", gen_number);
            for(int aux_x4=0; aux_x4<4; aux_x4++)
            {
                for(aux_y4=0; aux_y4<4; aux_y4++)
                {
                    for(aux_z4=0; aux_z4<4; aux_z4++)
                    {
                        printf("%d ", grid_odd[aux_x4][aux_y4][aux_z4]);
                    }
                    printf("\n");
                }
                printf("\n");
            }
            */

        }
        if(gen_number %2 == 0) // write to even
        {
            rules(number_of_cells, grid_even, grid_odd);
            /*
            printf("GENERATION %d\n", gen_number);
             for(int aux_x4=0; aux_x4<4; aux_x4++)
            {
                for(aux_y4=0; aux_y4<4; aux_y4++)
                {
                    for(aux_z4=0; aux_z4<4  ; aux_z4++)
                    {
                        printf("%d ", grid_even[aux_x4][aux_y4][aux_z4]);
                    }
                    printf("\n");
                }
                printf("\n");
            }
            */
        }

        //printf("1: %d current_gen: %d\n", count_species[0], gen_number);
        
        
        for(auxi=0; auxi < 9; auxi++)
        {
            if(count_species[auxi] > max_count[auxi])
            {
                max_count[auxi] = count_species[auxi];
                max_gen[auxi]=gen_number;
            }
        }
    }
    

    for(auxi=0; auxi < 9; auxi++)
    {
     
        printf("%ld %ld %ld \n", auxi+1, max_count[auxi], max_gen[auxi]);
    }

    liberarMatriz(number_of_cells);

    return 0;
}
