#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

int rank, size;
int my_coords[3];
int limite_inf_x, limite_inf_y , limite_inf_z ;
int limite_sup_x , limite_sup_y , limite_sup_z; 
float density;
int NUM_LINHAS;

unsigned int seed;
#define N_SPECIES 9

int *sub_divz_z;
int *sub_divz_y;
int *sub_divz_x;

char ***grid_even;
char ***grid_odd;

long count_species_local[10]={0,0,0,0,0,0,0,0,0,0};
long count_species[10]={0,0,0,0,0,0,0,0,0,0};
long max_count[10]={0,0,0,0,0,0,0,0,0,0};
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

void limites_x (){
    if (my_coords[2] == 0){
        limite_inf_x = 0; 
    }
        
    else{
        for (int i = 0; i <= my_coords[2] -1; ++ i){
            limite_inf_x = limite_inf_x + sub_divz_x [i];
        }
    }

    for (int i = 0; i<=my_coords[2]; ++i){
        limite_sup_x = limite_sup_x + sub_divz_x [i];
    }
}

void limites_y (){
    if (my_coords[1] == 0){
        limite_inf_y = 0; 
    }
        
    else{
        for (int i = 0; i <= my_coords[1] -1; ++ i){
            limite_inf_y = limite_inf_y + sub_divz_y [i];
        }
    }

    for (int i = 0; i<= my_coords[1]; ++i){
        limite_sup_y = limite_sup_y + sub_divz_y [i];
    }
}

void limites_z (){
    if (my_coords[0] == 0)
        limite_inf_z = 0;
    else{
        for (int i = 0; i <= my_coords[0] -1; ++ i){
            limite_inf_z = limite_inf_z + sub_divz_z [i];
        }
    }

    for (int i = 0; i<= my_coords[0]; ++i){
        limite_sup_z = limite_sup_z + sub_divz_z [i];
    }
}

void divide_number_parts(int number, int divide, int * sub_div) {

    int part_size, remainder;
    int start_index, end_index;
    int i;

    part_size = number / divide;
    remainder = number % divide;

    start_index = 0;

    for (i = 0; i < divide; i++) {
        end_index = start_index + part_size + (i < remainder ? 1 : 0);

        //printf("Part %d: ", i + 1);
        sub_div[i]=end_index-start_index;
        start_index = end_index;
    }
}

void My_MPI_Cart_Shift(MPI_Comm cart_comm, int pos_x, int pos_y,int pos_z, int dist_x, int dist_y,int dist_z, int *source, int*dest)
{
    MPI_Comm_rank(cart_comm, &rank);
    MPI_Cart_coords(cart_comm, rank, 3, my_coords);
    int my_coords1 = my_coords[pos_y];
    int my_coords2 = my_coords[pos_x];
    int my_coords3 = my_coords[pos_z];
    my_coords[pos_y]= my_coords1 + dist_y;
    my_coords[pos_x] = my_coords2 + dist_x;
    my_coords[pos_z] = my_coords3 + dist_z;
    MPI_Cart_rank(cart_comm, my_coords, dest);
    my_coords[pos_y]= my_coords1 -dist_y ;
    my_coords[pos_x]= my_coords2 -dist_x;
    my_coords[pos_z]= my_coords3 -dist_z;
    MPI_Cart_rank(cart_comm, my_coords, source);
}

void divide_em_tres (int *a_final, int *b_final, int *c_final, int size){
    int maior_linha = 0, maior_linha_prev = size;
    int a, b, c;
    int x = size;
    for (a = 1; a <= x; a++) {
        if (x % a == 0) {
            for (b = a; b <= x / a; b++) {
                if (x % (a * b) == 0) {
                    c = x / (a * b);
                    //printf("%d * %d * %d\n", a, b, c);
					
					if(maior_linha < a){
						maior_linha = a;
					}   
					
					if(maior_linha < b){
						maior_linha = b;
					}
					
					if(maior_linha < c){
						maior_linha = c;
					}
					
					//printf("maior_linha: %d\n", maior_linha);
                    // Atualizar os menores números encontrados até agora
                    if (maior_linha < maior_linha_prev) {
                      //  printf("entrou\n");
						*a_final = a;
                        *b_final = b;
                        *c_final = c;
						maior_linha_prev = maior_linha;
                    }
					maior_linha = 0;
                }
            }
        }
    }

}

void aloca_matrizes (int sub_x, int sub_y, int sub_z){
    grid_even = (char ***)malloc((sub_z+2) * sizeof(char **));
    grid_odd = (char ***)malloc((sub_z+2) * sizeof(char **));
    for (int i = 0; i < (sub_z+2); ++i) {
        grid_even[i] = (char **)malloc(((sub_y+2)) * sizeof(char *));
        grid_odd[i] = (char **)malloc(((sub_y+2)) * sizeof(char *));
        for (int j = 0; j < (sub_y+2); ++j) {
            grid_even[i][j] = (char *)malloc(((sub_x+2)) * sizeof(char));
            grid_odd[i][j] = (char *)malloc(((sub_x+2)) * sizeof(char)); 
        }
    }
}

void verifica_max (int gen_number,MPI_Comm MPI_COMM_WORLD){
    
    MPI_Reduce(count_species_local, count_species, sizeof (count_species), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    for (int x=1; x< 10; ++x){
        printf ("X=%d Cnt=%d\n", x, count_species[x]);
    }

    for(int x=1; x < 10; x++)
    {
        if(count_species[x] > max_count[x])
        {
            max_count[x] = count_species[x];
            max_gen[x]=gen_number;
        }
    }
}

void cria_primeira_grid (int NUM_LINHAS){
    int varrimento_x = 1;
    int varrimento_y = 1;
    int varrimento_z = 1;
    int valor_aux=0;

    int flag_y=0,flag_z=0;

    for (int init_z=0; init_z < NUM_LINHAS; init_z++){
        if (init_z >= limite_inf_z && init_z<limite_sup_z){
            flag_z = 1;
            ++varrimento_z;
        }
        else flag_z=0;

        for (int init_y=0; init_y < NUM_LINHAS; init_y++){
            if (init_y>=limite_inf_y && init_y<limite_sup_y){
                flag_y = 1;
                ++varrimento_y;
            }
            else flag_y=0;

            for (int init_x=0; init_x < NUM_LINHAS; init_x++){
                if(r4_uni() < density)
                {
                    // preenchimento initial do grid_even dependendo da seed
                    valor_aux = (int)(r4_uni() * N_SPECIES) + 1; // preenchimento initial do grid_even dependendo da see
                }else{
                  valor_aux = 0;
                }

                if (init_x>=limite_inf_x && init_x<limite_sup_x && flag_z == 1 && flag_y == 1 ){
                    grid_even[varrimento_z-1][varrimento_y-1][varrimento_x] = valor_aux;
                    count_species_local[valor_aux]++;
                    ++varrimento_x;
                }
            }
            //printf ("RANK :%d   Varrimento = %d\n",rank, varrimento_z);
            varrimento_x = 1;
        }
        varrimento_y = 1;
    }

    verifica_max (0);
}

void comunica_entre_processos (char ***data_send, int sub_x, int sub_y, int sub_z, MPI_Comm cart_comm){
    int aux_x, aux_y, aux_z; 
    int cima_rank, baixo_rank, esq_rank, dir_rank, frente_rank, tras_rank;
    int dir_cima_rank, esq_baixo_rank,dir_baixo_rank, esq_cima_rank, frente_cima_rank, tras_baixo_rank;
    int frente_baixo_rank, tras_cima_rank, dir_frente_rank, esq_tras_rank, dir_tras_rank, esq_frente_rank;
    int esq_cima_frente_rank, dir_baixo_tras_rank, dir_cima_frente_rank, esq_baixo_tras_rank;
    int esq_cima_tras_rank, dir_baixo_frente_rank, dir_cima_tras_rank, esq_baixo_frente_rank;

    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 0, 1, &esq_baixo_rank, &dir_cima_rank); // DIAG DIR CIMA/ ESQ BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 1, 0, &esq_tras_rank, &dir_frente_rank); // DIAG DIR CIMA/ ESQ BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, -1, 0, &esq_frente_rank, &dir_tras_rank); // DIAG DIR CIMA/ ESQ BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 0, -1, &esq_cima_rank, &dir_baixo_rank); // DIAG DIR BAIXO/ ESQ CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, -(NUM_LINHAS-1), 1, &tras_baixo_rank, &frente_cima_rank); // DIAG FRENTE CIMA / TRAS BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, -(NUM_LINHAS-1), -1, &tras_cima_rank, &frente_baixo_rank); // DIAG FRENTE CIMA / TRAS BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 0, 0, &esq_rank, &dir_rank); // FACE DIR/ESQ
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, 0, 1, &baixo_rank, &cima_rank); // FACE CIMA/BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, -1, 0, &frente_rank, &tras_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, -1, 1, 1, &dir_baixo_tras_rank, &esq_cima_frente_rank); // FACE TRAS/CIMA 
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 1, 1, &esq_baixo_tras_rank, &dir_cima_frente_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, -1, -1, 1, &dir_baixo_frente_rank, &esq_cima_tras_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, -1, 1, &esq_baixo_frente_rank, &dir_cima_tras_rank); // FACE TRAS/CIMA

    
    MPI_Sendrecv(&data_send[1][1][sub_x], 1, MPI_CHAR, dir_baixo_tras_rank, 0, &data_send[sub_z+1][sub_y+1][0], 1, MPI_CHAR, esq_cima_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    MPI_Sendrecv(&data_send[sub_z][sub_y][1], 1, MPI_CHAR, esq_cima_frente_rank, 0, &data_send[0][0][sub_x+1], 1, MPI_CHAR, dir_baixo_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[sub_z][sub_y][sub_x], 1, MPI_CHAR, dir_cima_frente_rank, 0, &data_send[0][0][0], 1, MPI_CHAR, esq_baixo_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[1][1][1], 1, MPI_CHAR, esq_baixo_tras_rank, 0, &data_send[sub_z+1][sub_y+1][sub_x+1], 1, MPI_CHAR, dir_cima_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    
    MPI_Sendrecv(&data_send[1][sub_y][sub_x], 1, MPI_CHAR, dir_baixo_frente_rank, 0, &data_send[sub_z+1][0][0], 1, MPI_CHAR, esq_cima_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    MPI_Sendrecv(&data_send[sub_z][1][1], 1, MPI_CHAR, esq_cima_tras_rank, 0, &data_send[0][sub_y+1][sub_x+1], 1, MPI_CHAR, dir_baixo_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[sub_z][1][sub_x], 1, MPI_CHAR, dir_cima_tras_rank, 0, &data_send[0][sub_y+1][0], 1, MPI_CHAR, esq_baixo_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[1][sub_y][1], 1, MPI_CHAR, esq_baixo_frente_rank, 0, &data_send[sub_z+1][0][sub_x+1], 1, MPI_CHAR, dir_cima_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    
    
    for( aux_z=0; aux_z < sub_z; aux_z++)
    {
        //FACE DIREITA DIAGS
        MPI_Sendrecv(&data_send[aux_z+1][1][1], 1, MPI_CHAR, esq_tras_rank, 0, &data_send[aux_z+1][sub_y+1][sub_x+1], 1, MPI_CHAR, dir_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
        MPI_Sendrecv(&data_send[aux_z+1][sub_y][1], 1, MPI_CHAR, esq_frente_rank, 0, &data_send[aux_z+1][0][sub_x+1], 1, MPI_CHAR, dir_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir baixo

        //FACE ESQUERDA DIAGS
        MPI_Sendrecv(&data_send[aux_z+1][1][sub_x], 1, MPI_CHAR, dir_tras_rank, 0, &data_send[aux_z+1][sub_y+1][0], 1, MPI_CHAR, esq_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[aux_z+1][sub_y][sub_x], 1, MPI_CHAR, dir_frente_rank, 0, &data_send[aux_z+1][0][0], 1, MPI_CHAR, esq_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        
        for (aux_y=0; aux_y<sub_y; aux_y++)
        {
            //FACE DIREITA 
            MPI_Sendrecv(&data_send[aux_z+1][aux_y+1][1], 1, MPI_CHAR, esq_rank, 0, &data_send[aux_z+1][aux_y+1][sub_x+1], 1, MPI_CHAR, dir_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
            
            //FACE ESQUERDA
            MPI_Sendrecv(&data_send[aux_z+1][aux_y+1][sub_x], 1, MPI_CHAR, dir_rank, 0, &data_send[aux_z+1][aux_y+1][0], 1, MPI_CHAR, esq_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir

        }
        for(aux_x=0; aux_x<sub_x; aux_x++)
        {
            //FACE FRENTE
            MPI_Sendrecv(&data_send[aux_z+1][1][aux_x+1], 1, MPI_CHAR, tras_rank, 0, &data_send[aux_z+1][sub_y+1][aux_x+1], 1, MPI_CHAR, frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir   
            
            //FACE TRAS
            MPI_Sendrecv(&data_send[aux_z+1][sub_y][aux_x+1], 1, MPI_CHAR, frente_rank, 0, &data_send[aux_z+1][0][aux_x+1], 1, MPI_CHAR, tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir   
        }
        
    }
    
    for (aux_y=0; aux_y<sub_y; aux_y++)
    { 
        //FACE CIMA DIAGS
        MPI_Sendrecv(&data_send[1][aux_y+1][sub_x],1, MPI_CHAR, dir_baixo_rank, 0, &data_send[sub_z+1][aux_y+1][0], 1, MPI_CHAR, esq_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[1][aux_y+1][1], 1, MPI_CHAR, esq_baixo_rank, 0, &data_send[sub_z+1][aux_y+1][sub_x+1], 1, MPI_CHAR, dir_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    
        //FACE BAIXO DIAGS
        MPI_Sendrecv(&data_send[sub_z][aux_y+1][sub_x], 1, MPI_CHAR, dir_cima_rank, 0, &data_send[0][aux_y+1][0], 1, MPI_CHAR, esq_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[sub_z][aux_y+1][1], 1, MPI_CHAR, esq_cima_rank, 0, &data_send[0][aux_y+1][sub_x+1], 1, MPI_CHAR, dir_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima

        for(aux_x=0; aux_x<sub_x; aux_x++)
        {
            //FACE CIMA
            MPI_Sendrecv(&data_send[1][aux_y+1][aux_x+1], 1, MPI_CHAR, baixo_rank, 0, &data_send[sub_z+1][aux_y+1][aux_x+1], 1, MPI_CHAR, cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir  
            //if(rank==4) printf("DATA SEND %d\n", data_send[1][aux_y+1][aux_x+1]);
            //FACE BAIXO
            MPI_Sendrecv(&data_send[sub_z][aux_y+1][aux_x+1], 1, MPI_CHAR, cima_rank, 0, &data_send[0][aux_y+1][aux_x+1], 1, MPI_CHAR, baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir                      
        }
    }
    
    for(aux_x=0; aux_x<sub_x; aux_x++)
    {
        //FACE FRENTE DIAGS
        MPI_Sendrecv(&data_send[sub_z][1][aux_x+1], 1, MPI_CHAR, tras_cima_rank, 0, &data_send[0][sub_y+1][aux_x+1], 1, MPI_CHAR, frente_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[1][1][aux_x+1], 1, MPI_CHAR, tras_baixo_rank, 0, &data_send[sub_z+1][sub_y+1][aux_x+1], 1, MPI_CHAR, frente_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima

        //FACE TRAS DIAGS
        MPI_Sendrecv(&data_send[sub_z][sub_y][aux_x+1], 1, MPI_CHAR, frente_cima_rank, 0, &data_send[0][0][aux_x+1], 1, MPI_CHAR, tras_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[1][sub_y][aux_x+1], 1, MPI_CHAR, frente_baixo_rank, 0, &data_send[sub_z+1][0][aux_x+1], 1, MPI_CHAR, tras_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    }
}

/************************************************************************************************
* Nome:death_rule
* funcao: verifica os vizinhos no caso da celula ter estado morta na geracao anterior
*
************************************************************************************************/
int death_rule(char *** grid, long aux_x, long aux_y, long aux_z)
{

    long search_x, search_y, search_z;
    long aux_search_y, aux_search_z;
    long cont_species_death[9]={0,0,0,0,0,0,0,0,0};
    int cont_rule=0;
    int max=0, max_pos=0, i;
    int x,y,z;
    
    for(search_x= aux_x-1, x=0; x < 3; x++, search_x++) 
    {
        for(search_y=aux_y-1, y=0; y < 3; y++, search_y++)
        {
            for(search_z=aux_z-1, z=0; z< 3;z++, search_z++)
            {
                if (grid [search_x][search_y][search_z] != 0){       
                    ++cont_rule;                
                    cont_species_death[grid[search_x][search_y][search_z]-1]++;
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
        for(i=1; i <9;i++ )
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
int life_rule (char *** grid, long aux_x, long aux_y, long aux_z){
    long search_x, search_y, search_z;
    long aux_search_y, aux_search_z;
    int cont_rule=-1;
    int x,y,z;
    
    // corre vizinhos e em caso de extremo verifica o extremo oposto
    for(search_x= aux_x-1, x=0; x < 3; x++, search_x++) 
    {
        for(search_y= aux_y-1, y=0; y < 3; y++, search_y++)
        {
            for(search_z= aux_z-1, z=0; z< 3;z++, search_z++)
            {
                //verifica se o vizinho está vivo
                if (grid [search_x][search_y][search_z] != 0){       
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
void rules(int sub_x ,int sub_y, int sub_z , char ***grid_new, char ***grid_old)
{
    long aux_x, aux_y, aux_z;

    //#pragma omp parallel private (aux_y, aux_z)
    //{
        //#pragma omp for schedule (dynamic)
        
        for(aux_x=1; aux_x<= sub_x; aux_x ++)
        {
            for(aux_y=1; aux_y<= sub_y; aux_y++)
            {
                for(aux_z=1; aux_z<= sub_z; aux_z++)
                {
                    if(grid_old[aux_x][aux_y][aux_z]==0) // morto 
                    { 
                        grid_new[aux_x][aux_y][aux_z]= death_rule(grid_old, aux_x, aux_y, aux_z);
                    }
                    else
                    {  
                        grid_new[aux_x][aux_y][aux_z]= life_rule(grid_old, aux_x, aux_y, aux_z);     
                    }
                    // se a celula esta viva nesta geracao, aumentamos o numero no array contador 
                    count_species_local[grid_new[aux_x][aux_y][aux_z]]++;
                }
            }
        }
   // }
}


int main(int argc, char *argv[]) {

    int number_of_gens;

    number_of_gens = atoi (argv[1]);
    NUM_LINHAS = atoi (argv[2]);
    density = atof (argv[3]);
    seed = atoi (argv[4]);

    init_r4uni(seed);

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int a_final, b_final, c_final;
    divide_em_tres (&a_final, &b_final, &c_final, size);

    sub_divz_z= (int *)malloc( a_final * sizeof(int)); 
    sub_divz_y= (int *)malloc( b_final * sizeof(int)); 
    sub_divz_x= (int *)malloc( c_final* sizeof(int)); 

    divide_number_parts(NUM_LINHAS,  a_final, sub_divz_z);
    divide_number_parts(NUM_LINHAS,  b_final, sub_divz_y);
    divide_number_parts(NUM_LINHAS,  c_final, sub_divz_x);

  
    //printf("SIZE %d\n", size);
    int count=0;
    //Cartesiano : 
    int dims[3] = { a_final, b_final, c_final};  // ISTO TEM DE VIR DOS INTEIROS QUE MULTIPLICAM O Nº PROCESSO
    //ACHO QUE TEM DE SER SEMPRE OS MESMOS 3 EIXOS  
    int periods[3] = {1, 1, 1};  // Enable wraparound
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart_comm);


    MPI_Cart_coords(cart_comm, rank, 3, my_coords);
    int sub_z = sub_divz_z[my_coords[0]];
    int sub_y = sub_divz_y[my_coords[1]];
    int sub_x = sub_divz_x[my_coords[2]];

    aloca_matrizes (sub_x, sub_y, sub_z);

    limites_x ();
    limites_y ();
    limites_z();

    cria_primeira_grid (NUM_LINHAS);
    comunica_entre_processos (grid_even, sub_x, sub_y, sub_z, cart_comm);
    
    
    for (int gen_number = 1; gen_number< number_of_gens; ++ gen_number){
       
        for (int auxi = 0; auxi < 10; ++auxi){
            count_species[auxi]=0;  
        }

        if (gen_number % 2 == 1){
            rules (sub_x, sub_y, sub_z, grid_odd, grid_even);
        }   
        else{
            rules (sub_x, sub_y, sub_z, grid_even, grid_odd);
        }

        verifica_max (gen_number, MPI_COMM_WORLD);  
    }

    if (rank == 0){
        for(int auxi=1; auxi < 10; auxi++)
        {
            printf("%d %ld %d \n", auxi, max_count[auxi], max_gen[auxi]);
        }
    }

    //free (sub_y, sub_z)
    
    MPI_Finalize();
    return 0; 
}
