#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int rank, size;
int my_coords[3];
#define TAMANHO_GRID 2

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
    

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int count=0;
    //Cartesiano : 
    int dims[3] = {TAMANHO_GRID, TAMANHO_GRID, TAMANHO_GRID+1};  // 2x2x2 grid
    int periods[3] = {1, 1, 1};  // Enable wraparound
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart_comm);

    
    MPI_Cart_coords(cart_comm, rank, 3, my_coords);

    int cima_rank, baixo_rank, esq_rank, dir_rank, frente_rank, tras_rank;
    int dir_cima_rank, esq_baixo_rank,dir_baixo_rank, esq_cima_rank, frente_cima_rank, tras_baixo_rank;
    int frente_baixo_rank, tras_cima_rank, dir_frente_rank, esq_tras_rank, dir_tras_rank, esq_frente_rank;
    int esq_cima_frente_rank, dir_baixo_tras_rank, dir_cima_frente_rank, esq_baixo_tras_rank;
    int esq_cima_tras_rank, dir_baixo_frente_rank, dir_cima_tras_rank, esq_baixo_frente_rank;


    int ***data_send = (int ***)malloc((TAMANHO_GRID+2) * sizeof(int **));
    for (int i = 0; i < (TAMANHO_GRID+2); ++i) {
        data_send[i] = (int **)malloc((TAMANHO_GRID+2) * sizeof(int *));
        for (int j = 0; j < (TAMANHO_GRID+2); ++j) {
            data_send[i][j] = (int *)malloc((TAMANHO_GRID+2) * sizeof(int));
            for (int k = 0; k < (TAMANHO_GRID+2); ++k) {
                if(k!=0 && i!=0 && j!= 0 && k!= (TAMANHO_GRID+1) && i!= (TAMANHO_GRID+1) && j!= (TAMANHO_GRID+1)) 
                {
                    data_send[i][j][k]=rank*1000 +count;
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
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, TAMANHO_GRID-1, 1, &tras_baixo_rank, &frente_cima_rank); // DIAG FRENTE CIMA / TRAS BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, TAMANHO_GRID-1, -1, &tras_cima_rank, &frente_baixo_rank); // DIAG FRENTE CIMA / TRAS BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 0, 0, &esq_rank, &dir_rank); // FACE DIR/ESQ
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, 0, 1, &baixo_rank, &cima_rank); // FACE CIMA/BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, 1, 0, &frente_rank, &tras_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, -1, -1, 1, &dir_baixo_tras_rank, &esq_cima_frente_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, -1, 1, &esq_baixo_tras_rank, &dir_cima_frente_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, -1, 1, 1, &dir_baixo_frente_rank, &esq_cima_tras_rank); // FACE TRAS/CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 1, 1, &esq_baixo_frente_rank, &dir_cima_tras_rank); // FACE TRAS/CIMA

    



    MPI_Sendrecv(&data_send[1][1][TAMANHO_GRID], 1, MPI_INT, dir_baixo_tras_rank, 0, &data_send[TAMANHO_GRID+1][TAMANHO_GRID+1][0], 1, MPI_INT, esq_cima_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    MPI_Sendrecv(&data_send[TAMANHO_GRID][TAMANHO_GRID][1], 1, MPI_INT, esq_cima_frente_rank, 0, &data_send[0][0][TAMANHO_GRID+1], 1, MPI_INT, dir_baixo_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[TAMANHO_GRID][TAMANHO_GRID][TAMANHO_GRID], 1, MPI_INT, dir_cima_frente_rank, 0, &data_send[0][0][0], 1, MPI_INT, esq_baixo_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[1][1][1], 1, MPI_INT, esq_baixo_tras_rank, 0, &data_send[TAMANHO_GRID+1][TAMANHO_GRID+1][TAMANHO_GRID+1], 1, MPI_INT, dir_cima_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    
    MPI_Sendrecv(&data_send[1][TAMANHO_GRID][TAMANHO_GRID], 1, MPI_INT, dir_baixo_frente_rank, 0, &data_send[TAMANHO_GRID+1][0][0], 1, MPI_INT, esq_cima_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    MPI_Sendrecv(&data_send[TAMANHO_GRID][1][1], 1, MPI_INT, esq_cima_tras_rank, 0, &data_send[0][TAMANHO_GRID+1][TAMANHO_GRID+1], 1, MPI_INT, dir_baixo_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[TAMANHO_GRID][1][TAMANHO_GRID], 1, MPI_INT, dir_cima_tras_rank, 0, &data_send[0][TAMANHO_GRID+1][0], 1, MPI_INT, esq_baixo_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[1][TAMANHO_GRID][1], 1, MPI_INT, esq_baixo_frente_rank, 0, &data_send[TAMANHO_GRID+1][0][TAMANHO_GRID+1], 1, MPI_INT, dir_cima_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);

    for( aux_z=0; aux_z < TAMANHO_GRID; aux_z++)
    {

        //FACE DIREITA DIAGS
        MPI_Sendrecv(&data_send[aux_z+1][1][1], 1, MPI_INT, esq_tras_rank, 0, &data_send[aux_z+1][TAMANHO_GRID+1][TAMANHO_GRID+1], 1, MPI_INT, dir_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
        MPI_Sendrecv(&data_send[aux_z+1][TAMANHO_GRID][1], 1, MPI_INT, esq_frente_rank, 0, &data_send[aux_z+1][0][TAMANHO_GRID+1], 1, MPI_INT, dir_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir baixo

        //FACE ESQUERDA DIAGS
        MPI_Sendrecv(&data_send[aux_z+1][1][TAMANHO_GRID], 1, MPI_INT, dir_tras_rank, 0, &data_send[aux_z+1][TAMANHO_GRID+1][0], 1, MPI_INT, esq_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[aux_z+1][TAMANHO_GRID][TAMANHO_GRID], 1, MPI_INT, dir_frente_rank, 0, &data_send[aux_z+1][0][0], 1, MPI_INT, esq_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo

        //FACE CIMA DIAGS
        MPI_Sendrecv(&data_send[1][aux_z+1][TAMANHO_GRID],1, MPI_INT, dir_baixo_rank, 0, &data_send[TAMANHO_GRID+1][aux_z+1][0], 1, MPI_INT, esq_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[1][aux_z+1][1], 1, MPI_INT, esq_baixo_rank, 0, &data_send[TAMANHO_GRID+1][aux_z+1][TAMANHO_GRID+1], 1, MPI_INT, dir_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima

        //FACE BAIXO DIAGS
        MPI_Sendrecv(&data_send[TAMANHO_GRID][aux_z+1][TAMANHO_GRID], 1, MPI_INT, dir_cima_rank, 0, &data_send[0][aux_z+1][0], 1, MPI_INT, esq_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[TAMANHO_GRID][aux_z+1][1], 1, MPI_INT, esq_cima_rank, 0, &data_send[0][aux_z+1][TAMANHO_GRID+1], 1, MPI_INT, dir_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima

        //FACE FRENTE DIAGS
        MPI_Sendrecv(&data_send[TAMANHO_GRID][1][aux_z+1], 1, MPI_INT, tras_cima_rank, 0, &data_send[0][TAMANHO_GRID+1][aux_z+1], 1, MPI_INT, frente_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[1][1][aux_z+1], 1, MPI_INT, tras_baixo_rank, 0, &data_send[TAMANHO_GRID+1][TAMANHO_GRID+1][aux_z+1], 1, MPI_INT, frente_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima

        //FACE TRAS DIAGS
        MPI_Sendrecv(&data_send[TAMANHO_GRID][TAMANHO_GRID][aux_z+1], 1, MPI_INT, frente_cima_rank, 0, &data_send[0][0][aux_z+1], 1, MPI_INT, tras_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[1][TAMANHO_GRID][aux_z+1], 1, MPI_INT, frente_baixo_rank, 0, &data_send[TAMANHO_GRID+1][0][aux_z+1], 1, MPI_INT, tras_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima

        

        for (aux_y=0; aux_y<TAMANHO_GRID; aux_y++)
        {
            //FACE DIREITA 
            MPI_Sendrecv(&data_send[aux_z+1][aux_y+1][1], 1, MPI_INT, esq_rank, 0, &data_send[aux_z+1][aux_y+1][TAMANHO_GRID+1], 1, MPI_INT, dir_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir

            //FACE ESQUERDA
            MPI_Sendrecv(&data_send[aux_z+1][aux_y+1][TAMANHO_GRID], 1, MPI_INT, dir_rank, 0, &data_send[aux_z+1][aux_y+1][0], 1, MPI_INT, esq_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
            
            //FACE CIMA
            MPI_Sendrecv(&data_send[1][aux_z+1][aux_y+1], 1, MPI_INT, baixo_rank, 0, &data_send[TAMANHO_GRID+1][aux_z+1][aux_y+1], 1, MPI_INT, cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir  
            
            //FACE BAIXO
            MPI_Sendrecv(&data_send[TAMANHO_GRID][aux_z+1][aux_y+1], 1, MPI_INT, cima_rank, 0, &data_send[0][aux_z+1][aux_y+1], 1, MPI_INT, baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir                      

            //FACE FRENTE
            MPI_Sendrecv(&data_send[aux_z+1][1][aux_y+1], 1, MPI_INT, tras_rank, 0, &data_send[aux_z+1][TAMANHO_GRID+1][aux_y+1], 1, MPI_INT, frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir   
            
            //FACE TRAS
            MPI_Sendrecv(&data_send[aux_z+1][TAMANHO_GRID][aux_y+1], 1, MPI_INT, frente_rank, 0, &data_send[aux_z+1][0][aux_y+1], 1, MPI_INT, tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir   

        }
        
    }

    if(rank==4)
    {
        
        for(aux_z=0;aux_z<(TAMANHO_GRID+2);aux_z++)
        {
            printf(" CAMADA\n");
            for(aux_y=0;aux_y<(TAMANHO_GRID+2);aux_y++)
            {
                for(aux_x=0;aux_x<(TAMANHO_GRID+2); aux_x++)
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
