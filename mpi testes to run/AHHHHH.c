#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int rank, size;
int my_coords[3];

void My_MPI_Cart_Shift(MPI_Comm cart_comm, int pos_x, int pos_y, int idk, int idk2, int *source, int*dest)
{
    MPI_Comm_rank(cart_comm, &rank);
    MPI_Cart_coords(cart_comm, rank, 3, my_coords);
    int my_coords1 = my_coords[idk];
    int my_coords2 = my_coords[idk2];
    my_coords[1]= my_coords1 + pos_x;
    my_coords[2] = my_coords2 + pos_y;
    MPI_Cart_rank(cart_comm, my_coords, dest);
    my_coords[1]= my_coords1 -pos_x ;
    my_coords[2]= my_coords2 -pos_y;
    MPI_Cart_rank(cart_comm, my_coords, source);
}

int main(int argc, char *argv[]) {
    

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Define Cartesian Topology
    int dims[3] = {3, 3, 3};  // 2x2x2 grid
    int periods[3] = {1, 1, 1};  // Enable wraparound
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart_comm);

    // Get Cartesian coordinates of current process
    
    //int my_coords_aux[3];
    MPI_Cart_coords(cart_comm, rank, 3, my_coords);
    // IMPORTANTE: COORD 0 - CIMA, COORD 1
    // Get neighbors
    int up_rank, down_rank, left_rank, right_rank, forward_rank, backward_rank, diag_rank, source_rank;
    MPI_Cart_shift(cart_comm, 0, 1, &up_rank, &down_rank);  
    MPI_Cart_shift(cart_comm, 1, 1, &left_rank, &right_rank);
    MPI_Cart_shift(cart_comm, 2, 1, &forward_rank, &backward_rank);
    MPI_Cart_shift(cart_comm, 1, 2, &source_rank, &diag_rank);
    //My_MPI_Cart_Shift(cart_comm, 1, 1, 2, 1, &source_rank, &diag_rank);

    //int diagonal_coords[2] = {(my_coords[0] + 1) % dims[0], (my_coords[1] + 1) % dims[1]};
    //int diagonal_rank_global;
    //MPI_Cart_rank(cart_comm, diagonal_coords, &diagonal_rank_global)



    // Send and receive data between neighbors
    //int data_send = rank * 10 + my_coords[0] * 100 + my_coords[1] * 1000 + my_coords[2] * 10000;  // Example data

    int ***data_send = (int ***)malloc(dims[0] * sizeof(int **));
    for (int i = 0; i < dims[0]; ++i) {
        data_send[i] = (int **)malloc(dims[1] * sizeof(int *));
        for (int j = 0; j < dims[1]; ++j) {
            data_send[i][j] = (int *)malloc(dims[2] * sizeof(int));
            for (int k = 0; k < dims[2]; ++k) {
                data_send[i][j][k]=rank;
                //data_send[i][j][k] = rank + my_coords[0] * 100 + my_coords[1] * 1000 + my_coords[2] * 10000;  // Example data
            }
        }
    }

    int ***data_recv_down = (int ***)malloc(dims[0] * sizeof(int **));
    for (int i = 0; i < dims[0]; ++i) {
        data_recv_down[i] = (int **)malloc(dims[1] * sizeof(int *));
        for (int j = 0; j < dims[1]; ++j) {
            data_recv_down[i][j] = (int *)malloc(dims[2] * sizeof(int));
        }
    }

    int aux;
    int data_recv_up, data_recv_left, data_recv_right, data_recv_forward, data_recv_backward;
    for( aux=0; aux < dims[2]; aux++)
    {
        //printf("rank: %d, SUPPOSED TO SEND %d\n",rank, data_send[0][aux]);
        MPI_Sendrecv(data_send[0][aux], dims[2], MPI_INT, source_rank, 0, data_recv_down[0][aux], dims[2], MPI_INT, diag_rank, 0, cart_comm, MPI_STATUS_IGNORE);
        //printf("Process %d, down: %d\n", rank, data_recv_down[0][aux]);
        
    }
    
    //MPI_Sendrecv(&data_send, 1, MPI_INT, down_rank, 0, &data_recv_up, 1, MPI_INT, up_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    //MPI_Sendrecv(&data_send, 1, MPI_INT, left_rank, 0, &data_recv_right, 1, MPI_INT, right_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    //MPI_Sendrecv(&data_send, 1, MPI_INT, right_rank, 0, &data_recv_left, 1, MPI_INT, left_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    //MPI_Sendrecv(&data_send, 1, MPI_INT, forward_rank, 0, &data_recv_backward, 1, MPI_INT, backward_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    //MPI_Sendrecv(&data_send, 1, MPI_INT, backward_rank, 0, &data_recv_forward, 1, MPI_INT, forward_rank, 0, cart_comm, MPI_STATUS_IGNORE);

    // Print received data
    if(rank==4)
    {
        for(aux=0;aux<dims[2];aux++)
        {
            printf("rank: %d, SUPPOSED TO SEND (FROM 0) %d\n",rank, data_send[0][0][aux]);
        }
    }

    if(rank==0)
    {
        for(aux=0;aux<dims[2];aux++)
        {
            printf("rank: %d, SUPPOSED TO RECEIVE (FROM 4) %d\n",rank, data_recv_down[0][0][aux]);
        }
    }

    MPI_Finalize();
    return 0;
}
