#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>

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
int *max_gen;
long *count_species;
long *count_species_new;

//long max_count[10]={0,0,0,0,0,0,0,0,0,0};


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
					
					if(maior_linha < a){
						maior_linha = a;
					}   
					
					if(maior_linha < b){
						maior_linha = b;
					}
					
					if(maior_linha < c){
						maior_linha = c;
					}
					
                    // Atualizar os menores números encontrados até agora
                    if (maior_linha < maior_linha_prev) {
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

char **alloc_2d_int(int rows, int cols) {
    /*
    char *data = (char *)calloc(rows*cols*sizeof(char));
    char **array= (char **)malloc(rows*sizeof(char*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);
    */
   char (* array)[cols]= malloc(sizeof((*array) *rows));
    return array;
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
            for(int k =0; k<(sub_x+2);k++)
            {
                grid_even[i][j][k]=0;
                grid_odd[i][j][k]=0;
            }
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
            varrimento_x = 1;
        }
        varrimento_y = 1;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(count_species_local, count_species, 10, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
if(rank==0)
{   
    for(int auxiii=0; auxiii < 10; auxiii++)
    {
        if(count_species[auxiii] > count_species_new[auxiii])
        {   
            count_species_new[auxiii] = count_species[auxiii];
            max_gen[auxiii]=0;
        }
    }    
}  


}

void freeMatrix2d(char ** matrix, int sub_y) {
    int i, j;

    for (i = 0; i < sub_y; i++) {
        free(matrix[i]);
    }
    free (matrix);

}

void comunica_entre_processos (char ***data_send, int sub_x, int sub_y, int sub_z,int coord_x, int coord_y, int coord_z, MPI_Comm cart_comm){
    
    //printf("SUB_X:%d    SUB_Y:%d    SUB_Z:%d\n", sub_x, sub_y, sub_z);

    char **face_dir_s, **face_dir_r, **face_esq_s, **face_esq_r, **face_cima_s, **face_cima_r, **face_baixo_s, **face_baixo_r;
    char **face_frente_s, **face_frente_r, **face_tras_s, **face_tras_r;
    char *diag_esq_tras_s, *diag_esq_tras_r, *diag_dir_tras_s, *diag_dir_tras_r, *diag_esq_frente_s, *diag_esq_frente_r, *diag_dir_frente_s, *diag_dir_frente_r;
    char *diag_esq_cima_s, *diag_esq_cima_r, *diag_dir_cima_s, *diag_dir_cima_r, *diag_esq_baixo_s, *diag_esq_baixo_r, *diag_dir_baixo_s, *diag_dir_baixo_r;
    char *diag_frente_baixo_s, *diag_frente_baixo_r, *diag_frente_cima_s, *diag_frente_cima_r, *diag_tras_baixo_s, *diag_tras_baixo_r, *diag_tras_cima_s, *diag_tras_cima_r;

    face_dir_s=alloc_2d_int(sub_z,sub_y);
    face_dir_r=alloc_2d_int(sub_z,sub_y);
    face_esq_s=alloc_2d_int(sub_z,sub_y);
    face_esq_r=alloc_2d_int(sub_z,sub_y);
    face_cima_s=alloc_2d_int(sub_y,sub_x);
    face_cima_r=alloc_2d_int(sub_y,sub_x);
    face_baixo_s=alloc_2d_int(sub_y,sub_x);
    face_baixo_r=alloc_2d_int(sub_y,sub_x);
    face_frente_s=alloc_2d_int(sub_z,sub_x);
    face_frente_r=alloc_2d_int(sub_z,sub_x);
    face_tras_s=alloc_2d_int(sub_z,sub_x);
    face_tras_r=alloc_2d_int(sub_z,sub_x);


    diag_esq_tras_r = (char *)malloc(sub_z*sizeof(char));
    diag_esq_tras_s = (char *)malloc(sub_z*sizeof(char));
    diag_dir_tras_r = (char *)malloc(sub_z*sizeof(char));
    diag_dir_tras_s = (char *)malloc(sub_z*sizeof(char));
    diag_esq_frente_r = (char *)malloc(sub_z*sizeof(char));
    diag_esq_frente_s = (char *)malloc(sub_z*sizeof(char));
    diag_dir_frente_r = (char *)malloc(sub_z*sizeof(char));
    diag_dir_frente_s = (char *)malloc(sub_z*sizeof(char));

    diag_esq_cima_r = (char *)malloc(sub_y*sizeof(char));
    diag_esq_cima_s = (char *)malloc(sub_y*sizeof(char));
    diag_dir_cima_r = (char *)malloc(sub_y*sizeof(char));
    diag_dir_cima_s = (char *)malloc(sub_y*sizeof(char));
    diag_esq_baixo_r = (char *)malloc(sub_y*sizeof(char));
    diag_esq_baixo_s = (char *)malloc(sub_y*sizeof(char));
    diag_dir_baixo_r = (char *)malloc(sub_y*sizeof(char));
    diag_dir_baixo_s = (char *)malloc(sub_y*sizeof(char));

    diag_frente_cima_r = (char *)malloc(sub_x*sizeof(char));
    diag_frente_cima_s = (char *)malloc(sub_x*sizeof(char));
    diag_tras_cima_r = (char *)malloc(sub_x*sizeof(char));
    diag_tras_cima_s = (char *)malloc(sub_x*sizeof(char));
    diag_tras_baixo_r = (char *)malloc(sub_x*sizeof(char));
    diag_tras_baixo_s = (char *)malloc(sub_x*sizeof(char));
    diag_frente_baixo_r = (char *)malloc(sub_x*sizeof(char));
    diag_frente_baixo_s = (char *)malloc(sub_x*sizeof(char));


    int aux_x, aux_y, aux_z; 
    int cima_rank, baixo_rank, esq_rank, dir_rank, frente_rank, tras_rank;
    int dir_cima_rank, esq_baixo_rank,dir_baixo_rank, esq_cima_rank, frente_cima_rank, tras_baixo_rank;
    int frente_baixo_rank, tras_cima_rank, dir_frente_rank, esq_tras_rank, dir_tras_rank, esq_frente_rank;
    int esq_cima_frente_rank, dir_baixo_tras_rank, dir_cima_frente_rank, esq_baixo_tras_rank;
    int esq_cima_tras_rank, dir_baixo_frente_rank, dir_cima_tras_rank, esq_baixo_frente_rank;

    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, coord_x+1, 0, coord_z+1, &esq_baixo_rank, &dir_cima_rank); // DIAG DIR CIMA/ ESQ BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, coord_x+1, -(coord_y-1), 0, &esq_tras_rank, &dir_frente_rank); // DIAG DIR CIMA/ ESQ BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, coord_x+1, -(coord_y+1), 0, &esq_frente_rank, &dir_tras_rank); // DIAG DIR CIMA/ ESQ BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, coord_x+1, 0, coord_z-1, &esq_cima_rank, &dir_baixo_rank); // DIAG DIR BAIXO/ ESQ CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, -(coord_y-1), coord_z+1, &tras_baixo_rank, &frente_cima_rank); // DIAG FRENTE CIMA / TRAS BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, -(coord_y-1), coord_z-1, &tras_cima_rank, &frente_baixo_rank); // DIAG FRENTE CIMA / TRAS BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, coord_x+1, 0, 0, &esq_rank, &dir_rank); // FACE DIR/ESQ
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, 0, coord_z+1, &baixo_rank, &cima_rank); // FACE CIMA/BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, -(coord_y-1), 0, &frente_rank, &tras_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, coord_x-1, -(coord_y-1), coord_z+1, &dir_baixo_tras_rank, &esq_cima_frente_rank); // FACE TRAS/CIMA 
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, coord_x+1, -(coord_y-1), coord_z+1, &esq_baixo_tras_rank, &dir_cima_frente_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, coord_x-1, -(coord_y+1), coord_z+1, &dir_baixo_frente_rank, &esq_cima_tras_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, coord_x+1, -(coord_y+1), coord_z+1, &esq_baixo_frente_rank, &dir_cima_tras_rank); // FACE TRAS/CIMA

    
    MPI_Sendrecv(&data_send[1][1][sub_x], 1, MPI_CHAR, dir_baixo_tras_rank, 0, &data_send[sub_z+1][sub_y+1][0], 1, MPI_CHAR, esq_cima_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    MPI_Sendrecv(&data_send[sub_z][sub_y][1], 1, MPI_CHAR, esq_cima_frente_rank, 0, &data_send[0][0][sub_x+1], 1, MPI_CHAR, dir_baixo_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[sub_z][sub_y][sub_x], 1, MPI_CHAR, dir_cima_frente_rank, 0, &data_send[0][0][0], 1, MPI_CHAR, esq_baixo_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[1][1][1], 1, MPI_CHAR, esq_baixo_tras_rank, 0, &data_send[sub_z+1][sub_y+1][sub_x+1], 1, MPI_CHAR, dir_cima_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    
    MPI_Sendrecv(&data_send[1][sub_y][sub_x], 1, MPI_CHAR, dir_baixo_frente_rank, 0, &data_send[sub_z+1][0][0], 1, MPI_CHAR, esq_cima_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    MPI_Sendrecv(&data_send[sub_z][1][1], 1, MPI_CHAR, esq_cima_tras_rank, 0, &data_send[0][sub_y+1][sub_x+1], 1, MPI_CHAR, dir_baixo_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[sub_z][1][sub_x], 1, MPI_CHAR, dir_cima_tras_rank, 0, &data_send[0][sub_y+1][0], 1, MPI_CHAR, esq_baixo_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[1][sub_y][1], 1, MPI_CHAR, esq_baixo_frente_rank, 0, &data_send[sub_z+1][0][sub_x+1], 1, MPI_CHAR, dir_cima_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    

    for(int k =0; k<sub_z;k++)
    {
        diag_dir_frente_s[k]=data_send[k+1][1][1];
        diag_dir_tras_s[k]=data_send[k+1][sub_y][1];

        diag_esq_frente_s[k]=data_send[k+1][1][sub_x];
        diag_esq_tras_s[k]=data_send[k+1][sub_y][sub_x];

        for(int i=0; i<sub_y; i++)
        {
            face_dir_s[k][i]=data_send[k+1][i+1][1];
            face_esq_s[k][i]=data_send[k+1][i+1][sub_x];
            //if (rank==1) printf("SEND %d\n", data_s[k][i]);
        }
        for(int j=0;j<sub_x;j++)
        {
            face_frente_s[k][j]=data_send[k+1][1][j+1];
            face_tras_s[k][j]=data_send[k+1][sub_y][j+1];
        }
        
    }

    for(int i=0; i<sub_y; i++)
    {
        diag_esq_cima_s[i]=data_send[1][i+1][sub_x];
        diag_dir_cima_s[i]=data_send[1][i+1][1];

        diag_esq_baixo_s[i]=data_send[sub_z][i+1][sub_x];
        diag_dir_baixo_s[i]=data_send[sub_z][i+1][1];

        for(int j=0;j<sub_x;j++)
        {
            face_cima_s[i][j]=data_send[1][i+1][j+1];
            face_baixo_s[i][j]=data_send[sub_z][i+1][j+1];
        }

    }

    for(int j=0;j<sub_x;j++)
    {
        diag_frente_baixo_s[j]=data_send[sub_z][1][j+1];
        diag_frente_cima_s[j]=data_send[1][1][j+1];

        diag_tras_baixo_s[j]=data_send[sub_z][sub_y][j+1];
        diag_tras_cima_s[j]=data_send[1][sub_y][j+1];
    }

    MPI_Sendrecv(&(diag_dir_frente_s[0]), sub_z, MPI_CHAR, esq_tras_rank, 0, &(diag_dir_frente_r[0]), sub_z, MPI_CHAR, dir_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_dir_tras_s[0]), sub_z, MPI_CHAR, esq_frente_rank, 0, &(diag_dir_tras_r[0]), sub_z, MPI_CHAR, dir_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_esq_frente_s[0]), sub_z, MPI_CHAR, dir_tras_rank, 0, &(diag_esq_frente_r[0]), sub_z, MPI_CHAR, esq_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_esq_tras_s[0]), sub_z, MPI_CHAR, dir_frente_rank, 0, &(diag_esq_tras_r[0]), sub_z, MPI_CHAR, esq_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    
    MPI_Sendrecv(&(diag_esq_cima_s[0]), sub_y, MPI_CHAR, dir_baixo_rank, 0, &(diag_esq_cima_r[0]), sub_y, MPI_CHAR, esq_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_dir_cima_s[0]), sub_y, MPI_CHAR, esq_baixo_rank, 0, &(diag_dir_cima_r[0]), sub_y, MPI_CHAR, dir_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_esq_baixo_s[0]), sub_y, MPI_CHAR, dir_cima_rank, 0, &(diag_esq_baixo_r[0]), sub_y, MPI_CHAR, esq_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_dir_baixo_s[0]), sub_y, MPI_CHAR, esq_cima_rank, 0, &(diag_dir_baixo_r[0]), sub_y, MPI_CHAR, dir_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    
    MPI_Sendrecv(&(diag_frente_baixo_s[0]), sub_x, MPI_CHAR, tras_cima_rank, 0, &(diag_frente_baixo_r[0]), sub_x, MPI_CHAR, frente_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_frente_cima_s[0]), sub_x, MPI_CHAR, tras_baixo_rank, 0, &(diag_frente_cima_r[0]), sub_x, MPI_CHAR, frente_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_tras_baixo_s[0]), sub_x, MPI_CHAR, frente_cima_rank, 0, &(diag_tras_baixo_r[0]), sub_x, MPI_CHAR, tras_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_tras_cima_s[0]), sub_x, MPI_CHAR, frente_baixo_rank, 0, &(diag_tras_cima_r[0]), sub_x, MPI_CHAR, tras_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
 
    MPI_Sendrecv(&(face_dir_s[0][0]), sub_z*sub_y, MPI_CHAR, esq_rank, 0, &(face_dir_r[0][0]), sub_z*sub_y, MPI_CHAR, dir_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(face_esq_s[0][0]), sub_z*sub_y, MPI_CHAR, dir_rank, 0, &(face_esq_r[0][0]), sub_z*sub_y, MPI_CHAR, esq_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(face_cima_s[0][0]), sub_y*sub_x, MPI_CHAR, baixo_rank, 0,&(face_cima_r[0][0]), sub_y*sub_x, MPI_CHAR, cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(face_baixo_s[0][0]), sub_y*sub_x, MPI_CHAR, cima_rank, 0, &(face_baixo_r[0][0]), sub_y*sub_x, MPI_CHAR, baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(face_frente_s[0][0]), sub_z*sub_x, MPI_CHAR, frente_rank, 0, &(face_frente_r[0][0]), sub_z*sub_x, MPI_CHAR, tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(face_tras_s[0][0]), sub_z*sub_x, MPI_CHAR, tras_rank, 0, &(face_tras_r[0][0]), sub_z*sub_x, MPI_CHAR, frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir

    for(int k =0; k<sub_z;k++)
    {

        data_send[k+1][sub_y+1][sub_x+1] = diag_dir_frente_r[k];
        data_send[k+1][0][sub_x+1] = diag_dir_tras_r[k];
        data_send[k+1][sub_y+1][0] = diag_esq_frente_r[k];
        data_send[k+1][0][0] = diag_esq_tras_r[k];

        for(int i=0; i<sub_y; i++)
        {
            data_send[k+1][i+1][sub_x+1]=face_dir_r[k][i];
            data_send[k+1][i+1][0]=face_esq_r[k][i];
            //if (rank==1) printf("SEND %d\n", data_s[k][i]);
        }
        for(int j=0;j<sub_x;j++)
        {
            data_send[k+1][sub_y+1][j+1]=face_frente_r[k][j];
            data_send[k+1][0][j+1]=face_tras_r[k][j];
        }
        
    }

    for(int i=0; i<sub_y; i++)
    {

        data_send[sub_z+1][i+1][0]=diag_esq_cima_r[i];
        data_send[sub_z+1][i+1][sub_x+1]=diag_dir_cima_r[i];
        data_send[0][i+1][0]=diag_esq_baixo_r[i];
        data_send[0][i+1][sub_x+1]=diag_dir_baixo_r[i];

        for(int j=0;j<sub_x;j++)
        {
            data_send[sub_z+1][i+1][j+1]=face_cima_r[i][j];
            data_send[0][i+1][j+1]=face_baixo_r[i][j];
        }

    }

    for(int j=0;j<sub_x;j++)
    {

        data_send[0][sub_y+1][j+1]=diag_frente_baixo_r[j];
        data_send[sub_z+1][sub_y+1][j+1]=diag_frente_cima_r[j];
        data_send[0][0][j+1]=diag_tras_baixo_r[j];
        data_send[sub_z+1][0][j+1]=diag_tras_cima_r[j];

    }


    free (face_dir_s[0]);
    free (face_dir_s);
    free (face_dir_r[0]);
    free (face_dir_r);
    free (face_esq_s[0]);
    free (face_esq_s);
    free (face_esq_r[0]);
    free (face_esq_r);
    free (face_cima_s[0]);
    free (face_cima_s);
    free (face_cima_r[0]);
    free (face_cima_r);
    free (face_baixo_s[0]);
    free (face_baixo_s);
    free (face_baixo_r[0]);
    free (face_baixo_r);
    free (face_frente_s[0]);
    free (face_frente_s);
    free (face_frente_r[0]);
    free (face_frente_r);
    free (face_tras_s[0]);
    free (face_tras_s);
    free (face_tras_r[0]);
    free (face_tras_r);
    

    free(diag_esq_tras_r);
    free(diag_esq_tras_s);
    free(diag_dir_tras_r);
    free(diag_dir_tras_s);
    free(diag_esq_frente_r);
    free(diag_esq_frente_s);
    free(diag_dir_frente_r);
    free(diag_dir_frente_s);

    free(diag_esq_cima_r);
    free(diag_esq_cima_s);
    free(diag_dir_cima_r);
    free(diag_dir_cima_s);
    free(diag_esq_baixo_r);
    free(diag_esq_baixo_s);
    free(diag_dir_baixo_r);
    free(diag_dir_baixo_s);

    free(diag_frente_cima_r);
    free(diag_frente_cima_s);
    free(diag_tras_cima_r);
    free(diag_tras_cima_s);
    free(diag_tras_baixo_r);
    free(diag_tras_baixo_s);
    free(diag_frente_baixo_r);
    free(diag_frente_baixo_s);
    
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
    
    for(search_z= aux_z-1, z=0; z <3; z++, search_z++) 
    {
        for(search_y=aux_y-1, y=0; y < 3; y++, search_y++)
        {
            for(search_x=aux_x-1, x=0; x< 3;x++, search_x++)
            {
                if (grid[search_z][search_y][search_x] != 0){       
                    ++cont_rule;                
                    cont_species_death[grid[search_z][search_y][search_x]-1]++;
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
        for(i=0; i <9; i++)
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
    for(search_z= aux_z-1, z=0; z < 3; z++, search_z++) 
    {
        for(search_y= aux_y-1, y=0; y < 3; y++, search_y++)
        {
            for(search_x= aux_x-1, x=0; x< 3; x++, search_x++)
            {
                //verifica se o vizinho está vivo
                if (grid [search_z][search_y][search_x] != 0){       
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
        return grid[aux_z][aux_y][aux_x];
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

    #pragma omp parallel private (aux_y, aux_x)
    {
        #pragma omp for reduction(+ : count_species_local) schedule (dynamic)
        
        for(aux_z=1; aux_z<= sub_z; aux_z ++)
        {   
            for(aux_y=1; aux_y<= sub_y; aux_y++)
            {
                for(aux_x=1; aux_x<= sub_x; aux_x++)
                {
                    if(grid_old[aux_z][aux_y][aux_x]==0) // morto 
                    { 
                        grid_new[aux_z][aux_y][aux_x]= death_rule(grid_old, aux_x, aux_y, aux_z);
                    }
                    else
                    {  
                        grid_new[aux_z][aux_y][aux_x]= life_rule(grid_old, aux_x, aux_y, aux_z);     
                    }
                    // se a celula esta viva nesta geracao, aumentamos o numero no array contador 
                    count_species_local[grid_new[aux_z][aux_y][aux_x]]++;
                }
            }
        }
    }
}

void freeMatrix(int sub_y, int sub_z) {
    int i, j;

    for (i = 0; i < sub_z; i++) {
        for (j = 0; j < sub_y; j++) {
            free(grid_even[i][j]);
            free(grid_odd[i][j]);
        }
        free(grid_even[i]);
        free(grid_odd[i]);
    }

    free(grid_even);
    free(grid_odd);
    free (max_gen);
    free (count_species);
    free (sub_divz_x);
    free (sub_divz_y);
    free (sub_divz_z);
    free (count_species_new);    
}


int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    double exec_time;

    max_gen = (int *)malloc( 10 * sizeof(int)); 
    count_species= (long *)malloc( 10 * sizeof(long)); 
    count_species_new= (long *)malloc( 10 * sizeof(long)); 
    
    for(int x=0; x<10; x++)
    {
        max_gen[x]=0;
        count_species[x]=0;
        count_species_new[x]=0;
    }

    int number_of_gens = 0;

    number_of_gens = atoi (argv[1]);
    NUM_LINHAS = atoi (argv[2]);
    density = atof (argv[3]);
    seed = atoi (argv[4]);

    init_r4uni(seed);

    int a_final, b_final, c_final;
    divide_em_tres (&a_final, &b_final, &c_final, size);
    //if(rank==0) printf("a :%d   b:%d    c:%d\n", a_final, b_final, c_final);

    sub_divz_z= (int *)malloc( a_final * sizeof(int)); 
    sub_divz_y= (int *)malloc( b_final * sizeof(int)); 
    sub_divz_x= (int *)malloc( c_final* sizeof(int)); 

    divide_number_parts(NUM_LINHAS,  a_final, sub_divz_z);
    divide_number_parts(NUM_LINHAS,  b_final, sub_divz_y);
    divide_number_parts(NUM_LINHAS,  c_final, sub_divz_x);

    


  
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
    if (rank==0)exec_time = -omp_get_wtime();

    comunica_entre_processos (grid_even, sub_x, sub_y, sub_z, c_final, b_final, a_final, cart_comm);

    for (int gen_number = 1; gen_number<= number_of_gens; ++gen_number){

        for (int auxi = 0; auxi < 10; ++auxi){
            count_species_local[auxi]=0;  
        }
        
        if (gen_number % 2 == 1){  
            rules (sub_x, sub_y, sub_z, grid_odd, grid_even);
            comunica_entre_processos (grid_odd, sub_x, sub_y, sub_z, c_final, b_final,a_final,  cart_comm);
        }   
        else{
            rules (sub_x, sub_y, sub_z, grid_even, grid_odd);
            comunica_entre_processos (grid_even, sub_x, sub_y, sub_z, c_final, b_final,a_final,  cart_comm);
        }
               
        MPI_Barrier(MPI_COMM_WORLD);
        
        for(int i=0; i<10;i++)
        {
             //printf("GEN:%d PROCESS: %d HAS LOCAL SPECIES %d with:%d\n", gen_number, rank,i, count_species_local[i]);
        }
 
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(count_species_local, count_species, 10, MPI_LONG, MPI_SUM,0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        if(rank==0)
        {
            //printf("GEN NUMBER:%d   \n", gen_number);
            for(int auxiii=1; auxiii < 10; auxiii++)
            {
                //printf("COUNT_SPECIES%d is %d    \n",auxiii, count_species[auxiii]);
                if(count_species[auxiii] > count_species_new[auxiii])
                {   
                    count_species_new[auxiii] = count_species[auxiii];
                    max_gen[auxiii]=gen_number;
                }
            }  
        }

        for (int auxi = 0; auxi < 10; ++auxi){
            count_species_local[auxi]=0;  
        }
        
    }

    if (rank==0){
        exec_time += omp_get_wtime();
        fprintf(stderr, "%.3fs\n", exec_time);
    } 

    if (rank == 0){
        for(int auxi=1; auxi < 10; auxi++)
        {
          printf("%d %ld %d \n", auxi, count_species_new[auxi], max_gen[auxi]);
        }
    }

    freeMatrix (sub_y, sub_z);

    

    
    
    MPI_Finalize();
    return 0; 
}
