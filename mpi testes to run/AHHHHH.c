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
    //IMPORATNTE: FACE DIREITA X=0,Y MUDA, Z MUDA
    //IMPORATNTE: FACE FRENTE X MUDA,Y=0, Z MUDA
    //IMPORTANTE: FACE ESQUERDA X=N, Y MUDA, Z MUDA
    //IMPORTANTE: FACE TRAS X MUDA, Y =N, Z MUDA
    //IMPORTANTE: FACE CIMA X MUDA, Y MUDA, Z=0
    //IMPORTANTE: FACE BAIXA X MUDA, Y MUDA, Z =N
    //IMPORTANTE: ARESTA CIMA DIREITA: X=0, Y MUDA, Z=0
    //IMPORTANTE: ARESTA CIMA FRENTE: X MUDA, Y=0, Z=0
    //IMPORTANTE: ARESTA CIMA ESQUERDA: X=N, Y=MUDA, Z=0
    //IMPORTANTE: ARESTA CIMA FRENTE: X MUDA, Y=N, Z=0
    //IMPORTANTE: ARESTA MEIO FRENTE_DIR: X=0, Y=0, Z MUDA
    //IMPORTANTE: ARESTA MEIO TRAS_DIR: X=0, Y=N, Z MUDA
    //IMPORTANTE: ARESTA MEIO FRENTE_ESQ: X=N, Y=0, Z MUDA
    //IMPORTANTE: ARESTA MEIO TRAS_ESQ: X=N, Y=N, Z MUDA
    //IMPORTANTE: ARESTA BAIXA DIREITA: X=0, Y MUDA, Z=N
    //IMPORTANTE: ARESTA BAIXA FRENTE: X MUDA, Y=0, Z=N
    //IMPORTANTE: ARESTA BAIXA ESQUERDA: X=N, Y=0, Z=N
    //IMPORTANTE: ARESTA BAIXA FRENTE: X MUDA, Y=N, Z=N
    // Get neighbors
    int cima_rank, baixo_rank, esq_rank, dir_rank, frente_rank, tras_rank;
    int dir_cima_rank, esq_baixo_rank,dir_baixo_rank, esq_cima_rank, frente_cima_rank, tras_baixo_rank;

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

    int ***data_send = (int ***)malloc(TAMANHO_GRID * sizeof(int **));
    for (int i = 0; i < TAMANHO_GRID; ++i) {
        data_send[i] = (int **)malloc(TAMANHO_GRID * sizeof(int *));
        for (int j = 0; j < TAMANHO_GRID; ++j) {
            data_send[i][j] = (int *)malloc(TAMANHO_GRID * sizeof(int));
            for (int k = 0; k < TAMANHO_GRID; ++k) {
                data_send[i][j][k]=rank*1000 +count;
                count++;
                //data_send[i][j][k] = rank + my_coords[0] * 100 + my_coords[1] * 1000 + my_coords[2] * 10000;  // Example data
            }
        }
    }

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


    int aux_x, aux_y, aux_z;
    int data_recv_up, data_recv_left, data_recv_right, data_recv_forward, data_recv_backward;
    //My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, 1, 1, &source_rank, &diag_rank); // frente cima
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 0, 1, &esq_baixo_rank, &dir_cima_rank); // DIAG DIR CIMA/ ESQ BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 0, -1, &esq_cima_rank, &dir_baixo_rank); // DIAG DIR BAIXO/ ESQ CIMA
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 0, TAMANHO_GRID-1, 1, &tras_baixo_rank, &frente_cima_rank); // DIAG FRENTE CIMA / TRAS BAIXO
    My_MPI_Cart_Shift(cart_comm, 2, 1, 0, 1, 0, 0, &esq_rank, &dir_rank); // FACE DIR/ESQ
    
    for( aux_z=0; aux_z < TAMANHO_GRID; aux_z++)
    {
        for (aux_y=0; aux_y<TAMANHO_GRID; aux_y++)
        {
            //FACE DIREITA INTEIRA
            MPI_Sendrecv(&data_send[aux_z][aux_y][0], dims[2], MPI_INT, esq_rank, 0, &data_recv_dir[aux_z+1][aux_y], dims[2], MPI_INT, dir_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir
            MPI_Sendrecv(&data_send[0][aux_y][0], dims[2], MPI_INT, esq_baixo_rank, 0, &data_recv_dir[TAMANHO_GRID+1][aux_y], dims[2], MPI_INT, dir_cima_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima
            MPI_Sendrecv(&data_send[TAMANHO_GRID-1][aux_y][0], dims[2], MPI_INT, esq_cima_rank, 0, &data_recv_dir[0][aux_y], dims[2], MPI_INT, dir_baixo_rank, 0, cart_comm, MPI_STATUS_IGNORE); // AR dir cima

            //FACE ESQUERDA
            MPI_Sendrecv(&data_send[aux_z][aux_y][TAMANHO_GRID-1], dims[2], MPI_INT, dir_rank, 0, &data_recv_esq[aux_z+1][aux_y], dims[2], MPI_INT, esq_rank, 0, cart_comm, MPI_STATUS_IGNORE); // face dir

            for(aux_x=0;aux_x<TAMANHO_GRID; aux_x++)
            {
                //
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
            for(aux_y=0;aux_y<TAMANHO_GRID;aux_y++)
            {
                printf("rank: %d, SUPPOSED TO RECEIVE FACE DIR (aux %d / %d) %d\n",rank,aux_z, aux_y, data_recv_dir[aux_z][aux_y]);
            }

            //printf("rank: %d, SUPPOSED TO RECEIVE FRENTE CIMA (aux %d) %d\n",rank,aux, data_recv_down[0][0][aux]);
            //printf("rank: %d, SUPPOSED TO RECEIVE TRAS BAIXO (aux %d) %d\n",rank,aux, data_recv_down[TAMANHO_GRID-1][TAMANHO_GRID-1][aux]);

            //printf("rank: %d, SUPPOSED TO RECEIVE DIR CIMA (aux %d) %d\n",rank,aux, data_recv_down[0][aux][0]);
            //printf("rank: %d, SUPPOSED TO RECEIVE ESQ BAIXO (aux %d) %d\n",rank,aux, data_recv_down[TAMANHO_GRID-1][aux][TAMANHO_GRID-1]);

        }

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
    }

    MPI_Finalize();
    return 0;
}
