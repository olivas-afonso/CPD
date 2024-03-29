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

//#define NUM_LINHAS 7


int *sub_divz_z;
int *sub_divz_y;
int *sub_divz_x;

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

int **alloc_2d_int(int rows, int cols) {
    int *data = (int *)malloc(rows*cols*sizeof(int));
    int **array= (int **)malloc(rows*sizeof(int*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

int main(int argc, char *argv[]) {

    int NUM_LINHAS;
    NUM_LINHAS= atoi (argv[1]);

    int varrimento_x = 1;
    int varrimento_y = 1;
    int varrimento_z = 1;
    int flag_y=0,flag_x=0;

    

    seed = 100;
    density=.4;
    init_r4uni(seed);

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


	int x = size;
    /*printf("Digite um inteiro x: ");
    scanf("%d", &x);
	*/ 
	
    //printf("Todas as combinações possíveis de três inteiros cujo produto é %d:\n", x);

    int maior_linha = 0, maior_linha_prev = x ;
    int a, b, c;
	int a_final, b_final, c_final;
	

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
						a_final = a;
                        b_final = b;
                        c_final = c;
						maior_linha_prev = maior_linha;
                    }
					maior_linha = 0;
                }
            }
        }
    }

    printf("A combinação com o maior menor número é: %d * %d * %d\n", a_final, b_final, c_final);

	
   
    // FALTA A DIVISAO DOS PROCESSADORES
    sub_divz_z= (int *)malloc( a_final * sizeof(int)); 
    /*
    for(int k=0; k<2; k++)
    {
        sub_divz_z[k]=aux_z_size[k];
    }
    */

    sub_divz_y= (int *)malloc( b_final * sizeof(int)); 

    /*
    for(int k=0; k<2; k++)
    {
        sub_divz_y[k]=aux_y_size[k];
    }
    */

    sub_divz_x= (int *)malloc( c_final* sizeof(int)); 

    /*
    for(int k=0; k<3; k++)
    {
        sub_divz_x[k]=aux_x_size[k];
    }
    */

	

   divide_number_parts(NUM_LINHAS,  a_final, sub_divz_z);
   divide_number_parts(NUM_LINHAS,  b_final, sub_divz_y);
   divide_number_parts(NUM_LINHAS,  c_final, sub_divz_x);

  
    //printf("SIZE %d\n", size);
    int count=0;
    //Cartesiano : 
    int dims[3] = { a_final, b_final, c_final};  // ISTO TEM DE VIR DOS INTEIROS QUE MULTIPLICAM O Nº PROCESSO
    //ACHO QUE TEM DE SER SEMPRE OS MESMOS 3 EIXOS ^^
    int periods[3] = {1, 1, 1};  // Enable wraparound
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart_comm);

    
    //int sub_divz_z[my_coords[0]] = NUM_LINHAS/SUB_DIV_Z;
    //int sub_divz_y[my_coords[1]] = SUB_DIV_Y;
    //int sub_divz_x[my_coords[2]] = SUB_DIV_X;

    MPI_Cart_coords(cart_comm, rank, 3, my_coords);
    int sub_z = sub_divz_z[my_coords[0]];
    int sub_y = sub_divz_y[my_coords[1]];
    int sub_x = sub_divz_x[my_coords[2]];

    int valor_aux=0;

    int ***data_send = (int ***)malloc((sub_z+2) * sizeof(int **));
    for (int i = 0; i < (sub_z+2); ++i) {
        data_send[i] = (int **)malloc(((sub_y+2)) * sizeof(int *));
        for (int j = 0; j < (sub_y+2); ++j) {
            data_send[i][j] = (int *)malloc(((sub_x+2)) * sizeof(int));
            for (int k = 0; k < (sub_x+2); ++k) {
                if((k!=0) && (i!=0) && (j!= 0) && (k!= (sub_x+1)) && (i!= (sub_z+1)) && (j!= (sub_y+1)) )
                {
                    data_send[i][j][k]=rank*1000; 
                    data_send[i][j][k] = data_send[i][j][k] + count;
                    count++;
                }
                else
                {
                   data_send[i][j][k]=0; 
                } 
            }
        }
    }

    int **send_x = (int **) malloc ((sub_z) * sizeof (int*));
    for(int i=0;i<sub_z;i++)
    {
        send_x[i] = (int *) malloc ((sub_y) * sizeof (int));
        for(int k =0; k<sub_y;k++)
        {
            send_x[i][k]=0;
        }
    }

    
    
    int **rcv_x = (int **) malloc ((sub_z) * sizeof (int*));
    for(int i=0;i<sub_z;i++)
    {
        rcv_x[i] = (int *) malloc ((sub_y) * sizeof (int));
        for(int k =0; k<sub_y;k++)
        {
            rcv_x[i][k]=0;
        }
    }

    
 
    int test[2][2];
    int **face_dir_s, **face_dir_r, **face_esq_s, **face_esq_r, **face_cima_s, **face_cima_r, **face_baixo_s, **face_baixo_r;
    int **face_frente_s, **face_frente_r, **face_tras_s, **face_tras_r;
    int *diag_esq_tras_s, *diag_esq_tras_r, *diag_dir_tras_s, *diag_dir_tras_r, *diag_esq_frente_s, *diag_esq_frente_r, *diag_dir_frente_s, *diag_dir_frente_r;
    int *diag_esq_cima_s, *diag_esq_cima_r, *diag_dir_cima_s, *diag_dir_cima_r, *diag_esq_baixo_s, *diag_esq_baixo_r, *diag_dir_baixo_s, *diag_dir_baixo_r;
    int *diag_frente_baixo_s, *diag_frente_baixo_r, *diag_frente_cima_s, *diag_frente_cima_r, *diag_tras_baixo_s, *diag_tras_baixo_r, *diag_tras_cima_s, *diag_tras_cima_r;
    
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


    diag_esq_tras_r = (int *)malloc(sub_z*sizeof(int));
    diag_esq_tras_s = (int *)malloc(sub_z*sizeof(int));
    diag_dir_tras_r = (int *)malloc(sub_z*sizeof(int));
    diag_dir_tras_s = (int *)malloc(sub_z*sizeof(int));
    diag_esq_frente_r = (int *)malloc(sub_z*sizeof(int));
    diag_esq_frente_s = (int *)malloc(sub_z*sizeof(int));
    diag_dir_frente_r = (int *)malloc(sub_z*sizeof(int));
    diag_dir_frente_s = (int *)malloc(sub_z*sizeof(int));

    diag_esq_cima_r = (int *)malloc(sub_y*sizeof(int));
    diag_esq_cima_s = (int *)malloc(sub_y*sizeof(int));
    diag_dir_cima_r = (int *)malloc(sub_y*sizeof(int));
    diag_dir_cima_s = (int *)malloc(sub_y*sizeof(int));
    diag_esq_baixo_r = (int *)malloc(sub_y*sizeof(int));
    diag_esq_baixo_s = (int *)malloc(sub_y*sizeof(int));
    diag_dir_baixo_r = (int *)malloc(sub_y*sizeof(int));
    diag_dir_baixo_s = (int *)malloc(sub_y*sizeof(int));

    diag_frente_cima_r = (int *)malloc(sub_x*sizeof(int));
    diag_frente_cima_s = (int *)malloc(sub_x*sizeof(int));
    diag_tras_cima_r = (int *)malloc(sub_x*sizeof(int));
    diag_tras_cima_s = (int *)malloc(sub_x*sizeof(int));
    diag_tras_baixo_r = (int *)malloc(sub_x*sizeof(int));
    diag_tras_baixo_s = (int *)malloc(sub_x*sizeof(int));
    diag_frente_baixo_r = (int *)malloc(sub_x*sizeof(int));
    diag_frente_baixo_s = (int *)malloc(sub_x*sizeof(int));

    int aux_x, aux_y, aux_z; 
    int cima_rank, baixo_rank, esq_rank, dir_rank, frente_rank, tras_rank;
    int dir_cima_rank, esq_baixo_rank,dir_baixo_rank, esq_cima_rank, frente_cima_rank, tras_baixo_rank;
    int frente_baixo_rank, tras_cima_rank, dir_frente_rank, esq_tras_rank, dir_tras_rank, esq_frente_rank;
    int esq_cima_frente_rank, dir_baixo_tras_rank, dir_cima_frente_rank, esq_baixo_tras_rank;
    int esq_cima_tras_rank, dir_baixo_frente_rank, dir_cima_tras_rank, esq_baixo_frente_rank;

    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, c_final+1, 0, a_final+1, &esq_baixo_rank, &dir_cima_rank); // DIAG DIR CIMA/ ESQ BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, c_final+1, b_final-1, 0, &esq_tras_rank, &dir_frente_rank); // DIAG DIR CIMA/ ESQ BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, c_final+1, b_final+1, 0, &esq_frente_rank, &dir_tras_rank); // DIAG DIR CIMA/ ESQ BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 0, a_final-1, &esq_cima_rank, &dir_baixo_rank); // DIAG DIR BAIXO/ ESQ CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, b_final-1, a_final+1, &tras_baixo_rank, &frente_cima_rank); // DIAG FRENTE CIMA / TRAS BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, b_final-1, a_final-1, &tras_cima_rank, &frente_baixo_rank); // DIAG FRENTE CIMA / TRAS BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, c_final+1, 0, 0, &esq_rank, &dir_rank); // FACE DIR/ESQ
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, 0, a_final+1, &baixo_rank, &cima_rank); // FACE CIMA/BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, b_final-1, 0, &frente_rank, &tras_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, c_final-1, b_final-1, a_final+1, &dir_baixo_tras_rank, &esq_cima_frente_rank); // FACE TRAS/CIMA 
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, c_final+1, b_final-1, a_final+1, &esq_baixo_tras_rank, &dir_cima_frente_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, c_final-1, b_final+1, a_final+1, &dir_baixo_frente_rank, &esq_cima_tras_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, c_final+1, b_final+1, a_final+1, &esq_baixo_frente_rank, &dir_cima_tras_rank); // FACE TRAS/CIMA

    
    MPI_Sendrecv(&data_send[1][1][sub_x], 1, MPI_INT, dir_baixo_tras_rank, 0, &data_send[sub_z+1][sub_y+1][0], 1, MPI_INT, esq_cima_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    MPI_Sendrecv(&data_send[sub_z][sub_y][1], 1, MPI_INT, esq_cima_frente_rank, 0, &data_send[0][0][sub_x+1], 1, MPI_INT, dir_baixo_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[sub_z][sub_y][sub_x], 1, MPI_INT, dir_cima_frente_rank, 0, &data_send[0][0][0], 1, MPI_INT, esq_baixo_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[1][1][1], 1, MPI_INT, esq_baixo_tras_rank, 0, &data_send[sub_z+1][sub_y+1][sub_x+1], 1, MPI_INT, dir_cima_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    
    MPI_Sendrecv(&data_send[1][sub_y][sub_x], 1, MPI_INT, dir_baixo_frente_rank, 0, &data_send[sub_z+1][0][0], 1, MPI_INT, esq_cima_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    MPI_Sendrecv(&data_send[sub_z][1][1], 1, MPI_INT, esq_cima_tras_rank, 0, &data_send[0][sub_y+1][sub_x+1], 1, MPI_INT, dir_baixo_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[sub_z][1][sub_x], 1, MPI_INT, dir_cima_tras_rank, 0, &data_send[0][sub_y+1][0], 1, MPI_INT, esq_baixo_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[1][sub_y][1], 1, MPI_INT, esq_baixo_frente_rank, 0, &data_send[sub_z+1][0][sub_x+1], 1, MPI_INT, dir_cima_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    

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

    MPI_Sendrecv(&(diag_dir_frente_s[0]), sub_z, MPI_INT, esq_tras_rank, 0, &(diag_dir_frente_r[0]), sub_z, MPI_INT, dir_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_dir_tras_s[0]), sub_z, MPI_INT, esq_frente_rank, 0, &(diag_dir_tras_r[0]), sub_z, MPI_INT, dir_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_esq_frente_s[0]), sub_z, MPI_INT, dir_tras_rank, 0, &(diag_esq_frente_r[0]), sub_z, MPI_INT, esq_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_esq_tras_s[0]), sub_z, MPI_INT, dir_frente_rank, 0, &(diag_esq_tras_r[0]), sub_z, MPI_INT, esq_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    
    MPI_Sendrecv(&(diag_esq_cima_s[0]), sub_y, MPI_INT, dir_baixo_rank, 0, &(diag_esq_cima_r[0]), sub_y, MPI_INT, esq_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_dir_cima_s[0]), sub_y, MPI_INT, esq_baixo_rank, 0, &(diag_dir_cima_r[0]), sub_y, MPI_INT, dir_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_esq_baixo_s[0]), sub_y, MPI_INT, dir_cima_rank, 0, &(diag_esq_baixo_r[0]), sub_y, MPI_INT, esq_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_dir_baixo_s[0]), sub_y, MPI_INT, esq_cima_rank, 0, &(diag_dir_baixo_r[0]), sub_y, MPI_INT, dir_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    
    MPI_Sendrecv(&(diag_frente_baixo_s[0]), sub_x, MPI_INT, tras_cima_rank, 0, &(diag_frente_baixo_r[0]), sub_x, MPI_INT, frente_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_frente_cima_s[0]), sub_x, MPI_INT, tras_baixo_rank, 0, &(diag_frente_cima_r[0]), sub_x, MPI_INT, frente_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_tras_baixo_s[0]), sub_x, MPI_INT, frente_cima_rank, 0, &(diag_tras_baixo_r[0]), sub_x, MPI_INT, tras_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(diag_tras_cima_s[0]), sub_x, MPI_INT, frente_baixo_rank, 0, &(diag_tras_cima_r[0]), sub_x, MPI_INT, tras_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
 
    MPI_Sendrecv(&(face_dir_s[0][0]), sub_z*sub_y, MPI_INT, esq_rank, 0, &(face_dir_r[0][0]), sub_z*sub_y, MPI_INT, dir_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(face_esq_s[0][0]), sub_z*sub_y, MPI_INT, dir_rank, 0, &(face_esq_r[0][0]), sub_z*sub_y, MPI_INT, esq_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(face_cima_s[0][0]), sub_y*sub_x, MPI_INT, cima_rank, 0,&(face_cima_r[0][0]), sub_y*sub_x, MPI_INT, baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(face_baixo_s[0][0]), sub_y*sub_x, MPI_INT, baixo_rank, 0, &(face_baixo_r[0][0]), sub_y*sub_x, MPI_INT, cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(face_frente_s[0][0]), sub_z*sub_x, MPI_INT, frente_rank, 0, &(face_frente_r[0][0]), sub_z*sub_x, MPI_INT, tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
    MPI_Sendrecv(&(face_tras_s[0][0]), sub_z*sub_x, MPI_INT, tras_rank, 0, &(face_tras_r[0][0]), sub_z*sub_x, MPI_INT, frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir

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
    
    for( aux_z=0; aux_z < sub_z; aux_z++)
    {
        /*
        //FACE DIREITA DIAGS
        MPI_Sendrecv(&data_send[aux_z+1][1][1], 1, MPI_INT, esq_tras_rank, 0, &data_send[aux_z+1][sub_y+1][sub_x+1], 1, MPI_INT, dir_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
        MPI_Sendrecv(&data_send[aux_z+1][sub_y][1], 1, MPI_INT, esq_frente_rank, 0, &data_send[aux_z+1][0][sub_x+1], 1, MPI_INT, dir_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir baixo

        //FACE ESQUERDA DIAGS
        MPI_Sendrecv(&data_send[aux_z+1][1][sub_x], 1, MPI_INT, dir_tras_rank, 0, &data_send[aux_z+1][sub_y+1][0], 1, MPI_INT, esq_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[aux_z+1][sub_y][sub_x], 1, MPI_INT, dir_frente_rank, 0, &data_send[aux_z+1][0][0], 1, MPI_INT, esq_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        */
        for (aux_y=0; aux_y<sub_y; aux_y++)
        {
            //FACE DIREITA 
            //MPI_Sendrecv(data_send[aux_z+1][aux_y+1][1], 1, MPI_INT, esq_rank, 0, &data_send[aux_z+1][aux_y+1][sub_x+1], 1, MPI_INT, dir_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
            
            //FACE ESQUERDA
            //MPI_Sendrecv(&data_send[aux_z+1][aux_y+1][sub_x], 1, MPI_INT, dir_rank, 0, &data_send[aux_z+1][aux_y+1][0], 1, MPI_INT, esq_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir

        }
        for(aux_x=0; aux_x<sub_x; aux_x++)
        {
            /*
            //FACE FRENTE
            MPI_Sendrecv(&data_send[aux_z+1][1][aux_x+1], 1, MPI_INT, tras_rank, 0, &data_send[aux_z+1][sub_y+1][aux_x+1], 1, MPI_INT, frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir   
            
            //FACE TRAS
            MPI_Sendrecv(&data_send[aux_z+1][sub_y][aux_x+1], 1, MPI_INT, frente_rank, 0, &data_send[aux_z+1][0][aux_x+1], 1, MPI_INT, tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir   
            */
        }
        
    }
    /*
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
    */
    //--------------------------------------DEBUG-----------------------------------------------
    if(rank==0)
    {
        
        printf("RANK: %d    SUB_Z: %d   SUB_Y: %d   SUB_X:  %d\n",rank, sub_z, sub_y, sub_x);
            //MPI_Cart_coords(cart_comm, rank, 3, my_coords);
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
