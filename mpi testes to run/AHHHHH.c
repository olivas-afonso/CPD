#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int rank, size;
int my_coords[3];
#define TAMANHO_GRID 3

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
    // Define Cartesian Topology
    int dims[3] = {3, 3, 3};  // 2x2x2 grid
    int periods[3] = {1, 1, 1};  // Enable wraparound
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart_comm);

    // Get Cartesian coordinates of current process
    
    //int my_coords_aux[3];
    MPI_Cart_coords(cart_comm, rank, 3, my_coords);
    // IMPORTANTE: COORD 0 - CAMADA, COORD 1 - COLUNA, COORD 2- LINHA
    //IMPORTANTE: NAO PRECISO DE VERTICES, MAS SIM DE ARESTAS + FACES.
    //IMPORTANTE: ACHO QUE AFINAL PRECISO DE VERTICES TAMBEM...
    //IMPORTANTE: POSSO DIMINUIR COMUNICACOES SE COPIAR VETORES ? SIM!!! VAMOS VER NO FUTURO SE ISTO VALE A PENA FAZER
    //IMPORTANTE: FALTA 2 DIAGONAIS EM CADA FACE.... COMO FAZER?? TENHO DE PENSAR!!!! AH!!!
    // Get neighbors
    int cima_rank, baixo_rank, esq_rank, dir_rank, frente_rank, tras_rank;
    int dir_cima_rank, esq_baixo_rank,dir_baixo_rank, esq_cima_rank, frente_cima_rank, tras_baixo_rank;
    int frente_baixo_rank, tras_cima_rank, dir_frente_rank, esq_tras_rank, dir_tras_rank, esq_frente_rank;
    int esq_cima_frente_rank, dir_baixo_tras_rank, dir_cima_frente_rank, esq_baixo_tras_rank;
    int esq_cima_tras_rank, dir_baixo_frente_rank, dir_cima_tras_rank, esq_baixo_frente_rank;

    //MPI_Cart_shift(cart_comm, 0, 1, &up_rank, &down_rank);  
    //MPI_Cart_shift(cart_comm, 1, 1, &left_rank, &right_rank);
    //MPI_Cart_shift(cart_comm, 2, 1, &forward_rank, &backward_rank);

    //MPI_Cart_shift(cart_comm, 2, 2, &source_rank, &diag_rank);
    //My_MPI_Cart_Shift(cart_comm, 1, 1, 2, 1, &source_rank, &diag_rank);
    //printf("ANTES\n");
    
    //My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 1, 1, &source_rank, &diag_rank);

    //int diagonal_coords[2] = {(my_coords[0] + 1) % dims[0], (my_coords[1] + 1) % dims[1]};
    //int diagonal_rank_global;
    //MPI_Cart_rank(cart_comm, diagonal_coords, &diagonal_rank_global)



    // Send and receive data between neighbors
    //int data_send = rank * 10 + my_coords[0] * 100 + my_coords[1] * 1000 + my_coords[2] * 10000;  // Example data

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
                //data_send[i][j][k] = rank + my_coords[0] * 100 + my_coords[1] * 1000 + my_coords[2] * 10000;  // Example data
            }
        }
    }
    /*
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==13)
    {
        for(int j=0; j< TAMANHO_GRID+2;j++)
        {
            for(int i=0; i< TAMANHO_GRID+2;i++)
            {
                for(int k=0; k< TAMANHO_GRID+2;k++)
                {
                    printf("DATA SEND j: %d i: %d  k: %d   : %d\n",j, i, k,  data_send[j][i][k]);
                }
            
            }
        }

    }
    */
    



    int ***data_recv_down = (int ***)malloc(TAMANHO_GRID * sizeof(int **));
    for (int i = 0; i < TAMANHO_GRID; ++i) {
        data_recv_down[i] = (int **)malloc(TAMANHO_GRID * sizeof(int *));
        for (int j = 0; j < TAMANHO_GRID; ++j) {
            data_recv_down[i][j] = (int *)malloc(TAMANHO_GRID * sizeof(int));
            for (int k = 0; k < TAMANHO_GRID; ++k) {
                data_recv_down[i][j][k]=0;
                
                //data_send[i][j][k] = rank + my_coords[0] * 100 + my_coords[1] * 1000 + my_coords[2] * 10000;  // Example data
            }
        }
    }

    int **data_recv_dir = (int **)malloc((TAMANHO_GRID+2) * sizeof(int *));
    for (int i = 0; i < (TAMANHO_GRID+2); ++i) {
        data_recv_dir[i] = (int *)malloc((TAMANHO_GRID) * sizeof(int));
        for (int j = 0; j < (TAMANHO_GRID+2); ++j) {
             data_recv_dir[i][j]=0;
        }
    }

    int **data_recv_esq = (int **)malloc((TAMANHO_GRID+2) * sizeof(int *));
    for (int i = 0; i < (TAMANHO_GRID+2); ++i) {
        data_recv_esq[i] = (int *)malloc((TAMANHO_GRID) * sizeof(int));
        for (int j = 0; j < (TAMANHO_GRID+2); ++j) {
             data_recv_esq[i][j]=0;
        }
    }

    int **data_recv_cima = (int **)malloc((TAMANHO_GRID+2) * sizeof(int *));
    for (int i = 0; i < (TAMANHO_GRID+2); ++i) {
        data_recv_cima[i] = (int *)malloc((TAMANHO_GRID) * sizeof(int));
        for (int j = 0; j < (TAMANHO_GRID+2); ++j) {
             data_recv_cima[i][j]=0;
        }
    }

    int **data_recv_baixo = (int **)malloc((TAMANHO_GRID+2) * sizeof(int *));
    for (int i = 0; i < (TAMANHO_GRID+2); ++i) {
        data_recv_baixo[i] = (int *)malloc((TAMANHO_GRID) * sizeof(int));
        for (int j = 0; j < (TAMANHO_GRID+2); ++j) {
             data_recv_baixo[i][j]=0;
        }
    }

    int **data_recv_frente = (int **)malloc((TAMANHO_GRID+2) * sizeof(int *));
    for (int i = 0; i < (TAMANHO_GRID+2); ++i) {
        data_recv_frente[i] = (int *)malloc((TAMANHO_GRID) * sizeof(int));
        for (int j = 0; j < (TAMANHO_GRID+2); ++j) {
             data_recv_frente[i][j]=0;
        }
    }

    int **data_recv_tras = (int **)malloc((TAMANHO_GRID+2) * sizeof(int *));
    for (int i = 0; i < (TAMANHO_GRID+2); ++i) {
        data_recv_tras[i] = (int *)malloc((TAMANHO_GRID) * sizeof(int));
        for (int j = 0; j < (TAMANHO_GRID+2); ++j) {
             data_recv_tras[i][j]=0;
        }
    }

    int vert_esq_cima_frente, vert_dir_baixo_tras, vert_dir_cima_frente, vert_esq_baixo_tras;
    int vert_esq_cima_tras, vert_dir_baixo_frente, vert_dir_cima_tras, vert_esq_baixo_frente;

    int aux_x, aux_y, aux_z;
    int data_recv_up, data_recv_left, data_recv_right, data_recv_forward, data_recv_backward;
    //My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, 1, 1, &source_rank, &diag_rank); // frente cima
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

    



    //VERT ESQ CIMA FRENTE
    if (rank==13) printf("START\n");
    MPI_Sendrecv(&data_send[1][1][TAMANHO_GRID], 1, MPI_INT, dir_baixo_tras_rank, 0, &vert_esq_cima_frente, 1, MPI_INT, esq_cima_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    MPI_Sendrecv(&data_send[TAMANHO_GRID][TAMANHO_GRID][1], 1, MPI_INT, esq_cima_frente_rank, 0, &vert_dir_baixo_tras, 1, MPI_INT, dir_baixo_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[TAMANHO_GRID][TAMANHO_GRID][TAMANHO_GRID], 1, MPI_INT, dir_cima_frente_rank, 0, &vert_esq_baixo_tras, 1, MPI_INT, esq_baixo_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[1][1][1], 1, MPI_INT, esq_baixo_tras_rank, 0, &vert_dir_cima_frente, 1, MPI_INT, dir_cima_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    
    MPI_Sendrecv(&data_send[1][TAMANHO_GRID][TAMANHO_GRID], 1, MPI_INT, dir_baixo_frente_rank, 0, &vert_esq_cima_tras, 1, MPI_INT, esq_cima_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
    MPI_Sendrecv(&data_send[TAMANHO_GRID][1][1], 1, MPI_INT, esq_cima_tras_rank, 0, &vert_dir_baixo_frente, 1, MPI_INT, dir_baixo_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[TAMANHO_GRID][1][TAMANHO_GRID], 1, MPI_INT, dir_cima_tras_rank, 0, &vert_esq_baixo_frente, 1, MPI_INT, esq_baixo_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send[1][TAMANHO_GRID][1], 1, MPI_INT, esq_baixo_frente_rank, 0, &vert_dir_cima_tras, 1, MPI_INT, dir_cima_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE);

    for( aux_z=0; aux_z < TAMANHO_GRID; aux_z++)
    {

        //FACE DIREITA DIAGS
        MPI_Sendrecv(&data_send[aux_z+1][1][1], 1, MPI_INT, esq_tras_rank, 0, &data_recv_dir[TAMANHO_GRID+1][aux_z], 1, MPI_INT, dir_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
        MPI_Sendrecv(&data_send[aux_z+1][TAMANHO_GRID][1], 1, MPI_INT, esq_frente_rank, 0, &data_recv_dir[0][aux_z], 1, MPI_INT, dir_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir baixo

        //FACE ESQUERDA DIAGS
        MPI_Sendrecv(&data_send[aux_z+1][1][TAMANHO_GRID], 1, MPI_INT, dir_tras_rank, 0, &data_recv_esq[TAMANHO_GRID+1][aux_z], 1, MPI_INT, esq_frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[aux_z+1][TAMANHO_GRID][TAMANHO_GRID], 1, MPI_INT, dir_frente_rank, 0, &data_recv_esq[0][aux_z], 1, MPI_INT, esq_tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo

        //FACE CIMA DIAGS
        MPI_Sendrecv(&data_send[1][aux_z+1][TAMANHO_GRID],1, MPI_INT, dir_baixo_rank, 0, &data_recv_cima[0][aux_z], 1, MPI_INT, esq_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[1][aux_z+1][1], 1, MPI_INT, esq_baixo_rank, 0, &data_recv_cima[TAMANHO_GRID+1][aux_z], 1, MPI_INT, dir_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima

        //FACE BAIXO DIAGS
        MPI_Sendrecv(&data_send[TAMANHO_GRID][aux_z+1][TAMANHO_GRID], 1, MPI_INT, dir_cima_rank, 0, &data_recv_baixo[0][aux_z], 1, MPI_INT, esq_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[TAMANHO_GRID][aux_z+1][1], 1, MPI_INT, esq_cima_rank, 0, &data_recv_baixo[TAMANHO_GRID+1][aux_z], 1, MPI_INT, dir_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima

        //FACE FRENTE DIAGS
        MPI_Sendrecv(&data_send[TAMANHO_GRID][1][aux_z+1], 1, MPI_INT, tras_cima_rank, 0, &data_recv_frente[0][aux_z], 1, MPI_INT, frente_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[1][1][aux_z+1], 1, MPI_INT, tras_baixo_rank, 0, &data_recv_frente[TAMANHO_GRID+1][aux_z], 1, MPI_INT, frente_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima

        //FACE TRAS DIAGS
        MPI_Sendrecv(&data_send[TAMANHO_GRID][TAMANHO_GRID][aux_z+1], 1, MPI_INT, frente_cima_rank, 0, &data_recv_tras[0][aux_z], 1, MPI_INT, tras_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo
        MPI_Sendrecv(&data_send[1][TAMANHO_GRID][aux_z+1], 1, MPI_INT, frente_baixo_rank, 0, &data_recv_tras[TAMANHO_GRID+1][aux_z], 1, MPI_INT, tras_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima

        

        for (aux_y=0; aux_y<TAMANHO_GRID; aux_y++)
        {
            //FACE DIREITA 
            MPI_Sendrecv(&data_send[aux_y+1][aux_z+1][1], 1, MPI_INT, esq_rank, 0, &data_send[aux_z+1][aux_y+1][TAMANHO_GRID+1], 1, MPI_INT, dir_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir

            //FACE ESQUERDA
            MPI_Sendrecv(&data_send[aux_y+1][aux_z+1][TAMANHO_GRID], 1, MPI_INT, dir_rank, 0, &data_recv_esq[aux_z+1][aux_y], 1, MPI_INT, esq_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
            
            //FACE CIMA
            MPI_Sendrecv(&data_send[1][aux_y+1][aux_z+1], 1, MPI_INT, baixo_rank, 0, &data_recv_cima[aux_z+1][aux_y], 1, MPI_INT, cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir  
            
            //FACE BAIXO
            MPI_Sendrecv(&data_send[TAMANHO_GRID][aux_y+1][aux_z+1], 1, MPI_INT, cima_rank, 0, &data_recv_baixo[aux_z+1][aux_y], 1, MPI_INT, baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir                      

            //FACE FRENTE
            MPI_Sendrecv(&data_send[aux_z+1][1][aux_y+1], 1, MPI_INT, tras_rank, 0, &data_recv_frente[aux_z+1][aux_y], 1, MPI_INT, frente_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir   
            
            //FACE TRAS
            MPI_Sendrecv(&data_send[aux_z+1][TAMANHO_GRID][aux_y+1], 1, MPI_INT, frente_rank, 0, &data_recv_tras[aux_z+1][aux_y], 1, MPI_INT, tras_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir   

            for(aux_x=0;aux_x<TAMANHO_GRID; aux_x++)
            {
                //FACE CIMA
                //MPI_Sendrecv(&data_send[0][aux_y][aux_x], dims[2], MPI_INT, baixo_rank, 0, &data_recv_cima[aux_y+1][aux_x], dims[2], MPI_INT, cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
                
            }
        }
        
        

        //MPI_Sendrecv(&data_send[0][0][aux], dims[2], MPI_INT, tras_baixo_rank, 0, &data_recv_down[0][0][aux], dims[2], MPI_INT, frente_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR fre cima
        //MPI_Sendrecv(&data_send[TAMANHO_GRID-1][TAMANHO_GRID-1][aux], dims[2], MPI_INT, frente_cima_rank, 0, &data_recv_down[TAMANHO_GRID-1][TAMANHO_GRID-1][aux], dims[2], MPI_INT, tras_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR tras baixo
        
        //MPI_Sendrecv(&data_send[0][aux][0], dims[2], MPI_INT, esq_baixo_rank, 0, &data_recv_down[0][aux][0], dims[2], MPI_INT, dir_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
        //MPI_Sendrecv(&data_send[TAMANHO_GRID-1][aux][TAMANHO_GRID-1], dims[2], MPI_INT, dir_cima_rank, 0, &data_recv_down[TAMANHO_GRID-1][aux][TAMANHO_GRID-1], dims[2], MPI_INT, esq_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR esq baixo


        //MPI_Sendrecv(&data_send[0][2][aux], dims[2], MPI_INT, source_rank, 0, &data_recv_down[0][2][aux], dims[2], MPI_INT, diag_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR cima fre
        /*
        MPI_Sendrecv(&data_send[dims[2]][aux][0], dims[2], MPI_INT, source_rank, 0, &data_recv_down[dims[2]][aux][0], dims[2], MPI_INT, diag_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR cima fre
        MPI_Sendrecv(&data_send[aux][0][0], dims[2], MPI_INT, source_rank, 0, &data_recv_down[aux][0][0], dims[2], MPI_INT, diag_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR cima fre
        MPI_Sendrecv(&data_send[aux][0][0], dims[2], MPI_INT, source_rank, 0, &data_recv_down[aux][0][0], dims[2], MPI_INT, diag_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR cima fre
        MPI_Sendrecv(&data_send[aux][0][0], dims[2], MPI_INT, source_rank, 0, &data_recv_down[aux][0][0], dims[2], MPI_INT, diag_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR cima fre
        MPI_Sendrecv(&data_send[aux][0][0], dims[2], MPI_INT, source_rank, 0, &data_recv_down[aux][0][0], dims[2], MPI_INT, diag_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR cima fre
        MPI_Sendrecv(&data_send[aux][0][0], dims[2], MPI_INT, source_rank, 0, &data_recv_down[aux][0][0], dims[2], MPI_INT, diag_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR cima fre
        */
        //printf("Process %d, down: %d\n", rank, data_recv_down[0][aux]);
        
    }
    
    //MPI_Sendrecv(&data_send, 1, MPI_INT, down_rank, 0, &data_recv_up, 1, MPI_INT, up_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    //MPI_Sendrecv(&data_send, 1, MPI_INT, left_rank, 0, &data_recv_right, 1, MPI_INT, right_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    //MPI_Sendrecv(&data_send, 1, MPI_INT, right_rank, 0, &data_recv_left, 1, MPI_INT, left_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    //MPI_Sendrecv(&data_send, 1, MPI_INT, forward_rank, 0, &data_recv_backward, 1, MPI_INT, backward_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    //MPI_Sendrecv(&data_send, 1, MPI_INT, backward_rank, 0, &data_recv_forward, 1, MPI_INT, forward_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);
    // Print received data


    if(rank==13)
    {
        for(aux_z=0;aux_z<(TAMANHO_GRID+2);aux_z++)
        {
            for(aux_y=0;aux_y<(TAMANHO_GRID+2);aux_y++)
            {
                printf("rank: %d, SUPPOSED TO RECEIVE FACE DIR (aux %d / %d) %d\n",rank,aux_z, aux_y, data_send[aux_z][aux_y][TAMANHO_GRID+1]);
            }

            //printf("rank: %d, SUPPOSED TO RECEIVE FRENTE CIMA (aux %d) %d\n",rank,aux, data_recv_down[0][0][aux]);
            //printf("rank: %d, SUPPOSED TO RECEIVE TRAS BAIXO (aux %d) %d\n",rank,aux, data_recv_down[TAMANHO_GRID-1][TAMANHO_GRID-1][aux]);

            //printf("rank: %d, SUPPOSED TO RECEIVE DIR CIMA (aux %d) %d\n",rank,aux, data_recv_down[0][aux][0]);
            //printf("rank: %d, SUPPOSED TO RECEIVE ESQ BAIXO (aux %d) %d\n",rank,aux, data_recv_down[TAMANHO_GRID-1][aux][TAMANHO_GRID-1]);
/*
        }
        printf("\n");
        for(aux_z=0;aux_z<(TAMANHO_GRID+2);aux_z++)
        {
            for(aux_y=0;aux_y<TAMANHO_GRID;aux_y++)
            {
                printf("rank: %d, SUPPOSED TO RECEIVE FACE ESQ (aux %d / %d) %d\n",rank,aux_z, aux_y, data_recv_esq[aux_z][aux_y]);
            }

            //printf("rank: %d, SUPPOSED TO RECEIVE FRENTE CIMA (aux %d) %d\n",rank,aux, data_recv_down[0][0][aux]);
            //printf("rank: %d, SUPPOSED TO RECEIVE TRAS BAIXO (aux %d) %d\n",rank,aux, data_recv_down[TAMANHO_GRID-1][TAMANHO_GRID-1][aux]);

            //printf("rank: %d, SUPPOSED TO RECEIVE DIR CIMA (aux %d) %d\n",rank,aux, data_recv_down[0][aux][0]);
            //printf("rank: %d, SUPPOSED TO RECEIVE ESQ BAIXO (aux %d) %d\n",rank,aux, data_recv_down[TAMANHO_GRID-1][aux][TAMANHO_GRID-1]);

        }
        printf("\n");
        for(aux_z=0;aux_z<(TAMANHO_GRID+2);aux_z++)
        {
            for(aux_y=0;aux_y<TAMANHO_GRID;aux_y++)
            {
                printf("rank: %d, SUPPOSED TO RECEIVE FACE CIMA (aux %d / %d) %d\n",rank,aux_z, aux_y, data_recv_cima[aux_z][aux_y]);
            }
        }

        printf("\n");
        for(aux_z=0;aux_z<(TAMANHO_GRID+2);aux_z++)
        {
            for(aux_y=0;aux_y<TAMANHO_GRID;aux_y++)
            {
                printf("rank: %d, SUPPOSED TO RECEIVE FACE BAIXO (aux %d / %d) %d\n",rank,aux_z, aux_y, data_recv_baixo[aux_z][aux_y]);
            }
        }

        printf("\n");
        for(aux_z=0;aux_z<(TAMANHO_GRID+2);aux_z++)
        {
            for(aux_y=0;aux_y<TAMANHO_GRID;aux_y++)
            {
                printf("rank: %d, SUPPOSED TO RECEIVE FACE FRENTE (aux %d / %d) %d\n",rank,aux_z, aux_y, data_recv_frente[aux_z][aux_y]);
            }
        }

        printf("\n");
        for(aux_z=0;aux_z<(TAMANHO_GRID+2);aux_z++)
        {
            for(aux_y=0;aux_y<TAMANHO_GRID;aux_y++)
            {
                printf("rank: %d, SUPPOSED TO RECEIVE FACE TRAS (aux %d / %d) %d\n",rank,aux_z, aux_y, data_recv_tras[aux_z][aux_y]);
            }
        }
        printf("\n");
        printf("rank: %d, SUPPOSED TO RECEIVE VERT ESQ CIMA FRENTE  %d\n",rank, vert_esq_cima_frente);

        printf("\n");
        printf("rank: %d, SUPPOSED TO RECEIVE VERT DIR BAIXO TRAS  %d\n",rank, vert_dir_baixo_tras);

        printf("\n");
        printf("rank: %d, SUPPOSED TO RECEIVE VERT ESQ BAIXO TRAS  %d\n",rank, vert_esq_baixo_tras);

        printf("\n");
        printf("rank: %d, SUPPOSED TO RECEIVE VERT DIR CIMA FRENTE  %d\n",rank, vert_dir_cima_frente);

        printf("\n");
        printf("rank: %d, SUPPOSED TO RECEIVE VERT ESQ CIMA TRAS  %d\n",rank, vert_esq_cima_tras);

        printf("\n");
        printf("rank: %d, SUPPOSED TO RECEIVE VERT DIR BAIXO FRENTE  %d\n",rank, vert_dir_baixo_frente);

        printf("\n");
        printf("rank: %d, SUPPOSED TO RECEIVE VERT ESQ BAIXO FRENTE  %d\n",rank, vert_esq_baixo_frente);

        printf("\n");
        printf("rank: %d, SUPPOSED TO RECEIVE VERT DIR CIMA TRAS  %d\n",rank, vert_dir_cima_tras);
*/}
}

   

    MPI_Finalize();
    return 0;
}
