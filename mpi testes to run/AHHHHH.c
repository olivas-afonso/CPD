#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int rank, size;
int my_coords[3];

//#define NUM_LINHAS 7


int *sub_divz_z;
int *sub_divz_y;
int *sub_divz_x;


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

int main(int argc, char *argv[]) {

    int NUM_LINHAS;
    NUM_LINHAS= atoi (argv[1]);
    //printf("OI\")
    
    // FALTA A DIVISAO DOS PROCESSADORES
    sub_divz_z= (int *)malloc(3 * sizeof(int)); 
    /*
    for(int k=0; k<2; k++)
    {
        sub_divz_z[k]=aux_z_size[k];
    }
    */

    sub_divz_y= (int *)malloc(3 * sizeof(int)); 

    /*
    for(int k=0; k<2; k++)
    {
        sub_divz_y[k]=aux_y_size[k];
    }
    */

    sub_divz_x= (int *)malloc(3* sizeof(int)); 

    /*
    for(int k=0; k<3; k++)
    {
        sub_divz_x[k]=aux_x_size[k];
    }
    */

   divide_number_parts(NUM_LINHAS, 2, sub_divz_z);
   divide_number_parts(NUM_LINHAS, 2, sub_divz_y);
   divide_number_parts(NUM_LINHAS, 3, sub_divz_x);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //printf("SIZE %d\n", size);
    int count=0;
    //Cartesiano : 
    int dims[3] = {3, 3, 3};  // ISTO TEM DE VIR DOS INTEIROS QUE MULTIPLICAM O Nº PROCESSO
    //ACHO QUE TEM DE SER SEMPRE OS MESMOS 3 EIXOS ^^
    int periods[3] = {1, 1, 1};  // Enable wraparound
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart_comm);

    

    int cima_rank, baixo_rank, esq_rank, dir_rank, frente_rank, tras_rank;
    int dir_cima_rank, esq_baixo_rank,dir_baixo_rank, esq_cima_rank, frente_cima_rank, tras_baixo_rank;
    int frente_baixo_rank, tras_cima_rank, dir_frente_rank, esq_tras_rank, dir_tras_rank, esq_frente_rank;
    int esq_cima_frente_rank, dir_baixo_tras_rank, dir_cima_frente_rank, esq_baixo_tras_rank;
    int esq_cima_tras_rank, dir_baixo_frente_rank, dir_cima_tras_rank, esq_baixo_frente_rank;

    //int sub_divz_z[my_coords[0]] = NUM_LINHAS/SUB_DIV_Z;
    //int sub_divz_y[my_coords[1]] = SUB_DIV_Y;
    //int sub_divz_x[my_coords[2]] = SUB_DIV_X;
    if(rank==0) printf("OI\n");
    MPI_Cart_coords(cart_comm, rank, 3, my_coords);
    int sub_z = sub_divz_z[my_coords[0]];
    int sub_y = sub_divz_y[my_coords[1]];
    int sub_x = sub_divz_x[my_coords[2]];

    /*
    if(rank==0)
    {
        printf("SUB_DIV_Z :%d   SUB_DIV_Z :%d \n",sub_divz_z[0],sub_divz_z[1]  );
        printf("SUB_DIV_Y :%d   SUB_DIV_Y :%d \n",sub_divz_y[0],sub_divz_y[1]  );
        printf("SUB_DIV_X :%d   SUB_DIV_X :%d   SUB_DIV_X :%d\n",sub_divz_x[0],sub_divz_x[1], sub_divz_x[2]  );
    }
    */

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

    
    if(rank==13)
    {
        //printf("RANK: %d    SUB_Z: %d   SUB_Y: %d   SUB_X:  %d\n",rank, sub_z, sub_y, sub_x);
            //MPI_Cart_coords(cart_comm, rank, 3, my_coords);
        for(aux_z=0;aux_z<((sub_z)+2);aux_z++)
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
