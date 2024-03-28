#include <stdio.h>

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int rank, size;
int my_coords[3];
int limite_inf_x, limite_inf_y , limite_inf_z ;
int limite_sup_x , limite_sup_y , limite_sup_z; 
float density;

unsigned int seed;
#define N_SPECIES 9

int *sub_divz_z;
int *sub_divz_y;
int *sub_divz_x;

char ***grid_even;
char ***grid_odd;

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

int main(int argc, char *argv[]) {
	
	printf("Comeca...\n");

    int NUM_LINHAS;
    int number_of_gens;

    number_of_gens = atoi (argv[1]);
    NUM_LINHAS = atoi (argv[2]);
    density = atof (argv[3]);
    seed = atoi (argv[4]);


    int varrimento_x = 1;
    int varrimento_y = 1;
    int varrimento_z = 1;
    int flag_y=0,flag_x=0;

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

    
    int cima_rank, baixo_rank, esq_rank, dir_rank, frente_rank, tras_rank;
    int dir_cima_rank, esq_baixo_rank,dir_baixo_rank, esq_cima_rank, frente_cima_rank, tras_baixo_rank;
    int frente_baixo_rank, tras_cima_rank, dir_frente_rank, esq_tras_rank, dir_tras_rank, esq_frente_rank;
    int esq_cima_frente_rank, dir_baixo_tras_rank, dir_cima_frente_rank, esq_baixo_tras_rank;
    int esq_cima_tras_rank, dir_baixo_frente_rank, dir_cima_tras_rank, esq_baixo_frente_rank;

    MPI_Cart_coords(cart_comm, rank, 3, my_coords);
    int sub_z = sub_divz_z[my_coords[0]];
    int sub_y = sub_divz_y[my_coords[1]];
    int sub_x = sub_divz_x[my_coords[2]];

    aloca_matrizes (sub_x, sub_y, sub_z);

    int valor_aux=0;

    limites_x ();
    limites_y ();
    limites_z();

    for (int init_x=0; init_x < NUM_LINHAS; init_x++){
    if (init_x >= limite_inf_z && init_x<limite_sup_z){
        flag_x = 1;
        ++varrimento_x;
    }
    else flag_x=0;
    
    for (int init_y=0; init_y < NUM_LINHAS; init_y++){
        if (init_y>=limite_inf_y && init_y<limite_sup_y){
            flag_y = 1;
            ++varrimento_y;
        }
        else flag_y=0;

        for (int init_z=0; init_z < NUM_LINHAS; init_z++){
            
             if(r4_uni() < density)
                    {
                        // preenchimento initial do grid_even dependendo da seed
                        valor_aux = (int)(r4_uni() * N_SPECIES) + 1; // preenchimento initial do grid_even dependendo da see
                    }else{
                        valor_aux = 0;
                    }

            if (init_z>=limite_inf_x && init_z<limite_sup_x && flag_x == 1 && flag_y == 1 ){

                grid_even[varrimento_x-1][varrimento_y-1][varrimento_z] = valor_aux;
                  //printf("VALORES A ENTRAR %d, pos_x = %d, pos_y = %d, pos_z = %d \n", data_send[varrimento_x-1][varrimento_y-1][varrimento_z], varrimento_x-1, varrimento_y-1, varrimento_z);
                 ++varrimento_z;
            }
        }
        //printf ("RANK :%d   Varrimento = %d\n",rank, varrimento_z);
        varrimento_z = 1;
    }
    
    varrimento_y = 1;
}

    
 
    

    int vert_esq_cima_frente, vert_dir_baixo_tras, vert_dir_cima_frente, vert_esq_baixo_tras;
    int vert_esq_cima_tras, vert_dir_baixo_frente, vert_dir_cima_tras, vert_esq_baixo_frente;

    int aux_x, aux_y, aux_z;
    int data_recv_up, data_recv_left, data_recv_right, data_recv_forward, data_recv_backward;
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 0, 1, &esq_baixo_rank, &dir_cima_rank); // DIAG DIR CIMA/ ESQ BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, -1, 0, &esq_tras_rank, &dir_frente_rank); // DIAG DIR CIMA/ ESQ BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 1, 0, &esq_frente_rank, &dir_tras_rank); // DIAG DIR CIMA/ ESQ BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 0, -1, &esq_cima_rank, &dir_baixo_rank); // DIAG DIR BAIXO/ ESQ CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, NUM_LINHAS-1, 1, &tras_baixo_rank, &frente_cima_rank); // DIAG FRENTE CIMA / TRAS BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, NUM_LINHAS-1, -1, &tras_cima_rank, &frente_baixo_rank); // DIAG FRENTE CIMA / TRAS BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 0, 0, &esq_rank, &dir_rank); // FACE DIR/ESQ
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, 0, 1, &baixo_rank, &cima_rank); // FACE CIMA/BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, 1, 0, &frente_rank, &tras_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, -1, -1, 1, &dir_baixo_tras_rank, &esq_cima_frente_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, -1, 1, &esq_baixo_tras_rank, &dir_cima_frente_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, -1, 1, 1, &dir_baixo_frente_rank, &esq_cima_tras_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 1, 1, &esq_baixo_frente_rank, &dir_cima_tras_rank); // FACE TRAS/CIMA

    


    
    MPI_Sendrecv(&data_send[1][1][sub_x], 1, MPI_INT, dir_baixo_tras_rank, 0, &data_send[sub_z+1][sub_y+1][0], 1, MPI_INT, esq_cima_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    MPI_Sendrecv(&data_send[sub_z][sub_y][1], 1, MPI_INT, esq_cima_frente_rank, 0, &data_send[0][0][sub_x+1], 1, MPI_INT, dir_baixo_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[sub_z][sub_y][sub_x], 1, MPI_INT, dir_cima_frente_rank, 0, &data_send[0][0][0], 1, MPI_INT, esq_baixo_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[1][1][1], 1, MPI_INT, esq_baixo_tras_rank, 0, &data_send[sub_z+1][sub_y+1][sub_x+1], 1, MPI_INT, dir_cima_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    
    MPI_Sendrecv(&data_send[1][sub_y][sub_x], 1, MPI_INT, dir_baixo_frente_rank, 0, &data_send[sub_z+1][0][0], 1, MPI_INT, esq_cima_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    MPI_Sendrecv(&data_send[sub_z][1][1], 1, MPI_INT, esq_cima_tras_rank, 0, &data_send[0][sub_y+1][sub_x+1], 1, MPI_INT, dir_baixo_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[sub_z][1][sub_x], 1, MPI_INT, dir_cima_tras_rank, 0, &data_send[0][sub_y+1][0], 1, MPI_INT, esq_baixo_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[1][sub_y][1], 1, MPI_INT, esq_baixo_frente_rank, 0, &data_send[sub_z+1][0][sub_x+1], 1, MPI_INT, dir_cima_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    
    
    for( aux_z=0; aux_z < sub_z; aux_z++)
    {
        
        //FACE DIREITA DIAGS
        MPI_Sendrecv(&data_send[aux_z+1][1][1], 1, MPI_INT, esq_tras_rank, 0, &data_send[aux_z+1][sub_y+1][sub_x+1], 1, MPI_INT, dir_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
        MPI_Sendrecv(&data_send[aux_z+1][sub_y][1], 1, MPI_INT, esq_frente_rank, 0, &data_send[aux_z+1][0][sub_x+1], 1, MPI_INT, dir_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir baixo

        //FACE ESQUERDA DIAGS
        MPI_Sendrecv(&data_send[aux_z+1][1][sub_x], 1, MPI_INT, dir_tras_rank, 0, &data_send[aux_z+1][sub_y+1][0], 1, MPI_INT, esq_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[aux_z+1][sub_y][sub_x], 1, MPI_INT, dir_frente_rank, 0, &data_send[aux_z+1][0][0], 1, MPI_INT, esq_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        
        
        for (aux_y=0; aux_y<sub_y; aux_y++)
        {
            //FACE DIREITA 
            MPI_Sendrecv(&data_send[aux_z+1][aux_y+1][1], 1, MPI_INT, esq_rank, 0, &data_send[aux_z+1][aux_y+1][sub_x+1], 1, MPI_INT, dir_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
            
            //FACE ESQUERDA
            MPI_Sendrecv(&data_send[aux_z+1][aux_y+1][sub_x], 1, MPI_INT, dir_rank, 0, &data_send[aux_z+1][aux_y+1][0], 1, MPI_INT, esq_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir

        }
        for(aux_x=0; aux_x<sub_x; aux_x++)
        {

            //FACE FRENTE
            MPI_Sendrecv(&data_send[aux_z+1][1][aux_x+1], 1, MPI_INT, tras_rank, 0, &data_send[aux_z+1][sub_y+1][aux_x+1], 1, MPI_INT, frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir   
            
            //FACE TRAS
            MPI_Sendrecv(&data_send[aux_z+1][sub_y][aux_x+1], 1, MPI_INT, frente_rank, 0, &data_send[aux_z+1][0][aux_x+1], 1, MPI_INT, tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir   
        }
        
    }
    

    for (aux_y=0; aux_y<sub_y; aux_y++)
    {
        
        //FACE CIMA DIAGS
        MPI_Sendrecv(&data_send[1][aux_y+1][sub_x],1, MPI_INT, dir_baixo_rank, 0, &data_send[sub_z+1][aux_y+1][0], 1, MPI_INT, esq_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[1][aux_y+1][1], 1, MPI_INT, esq_baixo_rank, 0, &data_send[sub_z+1][aux_y+1][sub_x+1], 1, MPI_INT, dir_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    
        //FACE BAIXO DIAGS
        MPI_Sendrecv(&data_send[sub_z][aux_y+1][sub_x], 1, MPI_INT, dir_cima_rank, 0, &data_send[0][aux_y+1][0], 1, MPI_INT, esq_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[sub_z][aux_y+1][1], 1, MPI_INT, esq_cima_rank, 0, &data_send[0][aux_y+1][sub_x+1], 1, MPI_INT, dir_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima

        for(aux_x=0; aux_x<sub_x; aux_x++)
        {
            //FACE CIMA
            MPI_Sendrecv(&data_send[1][aux_y+1][aux_x+1], 1, MPI_INT, baixo_rank, 0, &data_send[sub_z+1][aux_y+1][aux_x+1], 1, MPI_INT, cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir  
            //if(rank==4) printf("DATA SEND %d\n", data_send[1][aux_y+1][aux_x+1]);
            //FACE BAIXO
            MPI_Sendrecv(&data_send[sub_z][aux_y+1][aux_x+1], 1, MPI_INT, cima_rank, 0, &data_send[0][aux_y+1][aux_x+1], 1, MPI_INT, baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir                      
        }
    }
    
    for(aux_x=0; aux_x<sub_x; aux_x++)
    {
        //FACE FRENTE DIAGS
        MPI_Sendrecv(&data_send[sub_z][1][aux_x+1], 1, MPI_INT, tras_cima_rank, 0, &data_send[0][sub_y+1][aux_x+1], 1, MPI_INT, frente_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[1][1][aux_x+1], 1, MPI_INT, tras_baixo_rank, 0, &data_send[sub_z+1][sub_y+1][aux_x+1], 1, MPI_INT, frente_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima

        //FACE TRAS DIAGS
        MPI_Sendrecv(&data_send[sub_z][sub_y][aux_x+1], 1, MPI_INT, frente_cima_rank, 0, &data_send[0][0][aux_x+1], 1, MPI_INT, tras_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[1][sub_y][aux_x+1], 1, MPI_INT, frente_baixo_rank, 0, &data_send[sub_z+1][0][aux_x+1], 1, MPI_INT, tras_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    }

    //--------------------------------------DEBUG-----------------------------------------------
    if(rank==0)
    {
        //printf("RANK: %d    SUB_Z: %d   SUB_Y: %d   SUB_X:  %d\n",rank, sub_z, sub_y, sub_x);
            //MPI_Cart_coords(cart_comm, rank, 3, my_coords)
            printf("RANK :%d\n", rank);
        for(aux_z=0;aux_z<(sub_z+2);aux_z++)
        {
            printf("CAMADA %d\n", aux_z);
            for(aux_y=0;aux_y<(sub_y+2);aux_y++)
            {
                for(aux_x=0;aux_x<(sub_x+2); aux_x++)
                {
                     printf("%d ",data_send[aux_z][aux_y][aux_x]);
                }
                printf("\n");       
            }

        }
    }
    


    

   
    


    MPI_Finalize();
    return 0;
}
