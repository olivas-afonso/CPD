#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Define Cartesian Topology
    int dims[3] = {2, 2, 2};  // 2x2x2 grid
    int periods[3] = {1, 1, 1};  // Enable wraparound
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &cart_comm);

    // Get Cartesian coordinates of current process
    int my_coords[3];
    MPI_Cart_coords(cart_comm, rank, 3, my_coords);

    // Allocate memory for data_send
    int ***data_send = (int ***)malloc(dims[0] * sizeof(int **));
    for (int i = 0; i < dims[0]; ++i) {
        data_send[i] = (int **)malloc(dims[1] * sizeof(int *));
        for (int j = 0; j < dims[1]; ++j) {
            data_send[i][j] = (int *)malloc(dims[2] * sizeof(int));
            for (int k = 0; k < dims[2]; ++k) {
                data_send[i][j][k] = rank * 100 + i * 1000 + j * 10000 + k * 100000;  // Example data
            }
        }
    }

    // Get neighbors
    int up_rank, down_rank, left_rank, right_rank, forward_rank, backward_rank;
    MPI_Cart_shift(cart_comm, 0, 1, &up_rank, &down_rank);
    MPI_Cart_shift(cart_comm, 1, 1, &left_rank, &right_rank);
    MPI_Cart_shift(cart_comm, 2, 1, &forward_rank, &backward_rank);

    // Send and receive data between neighbors
    int ***data_recv_up = (int ***)malloc(dims[0] * sizeof(int **));
    int ***data_recv_down = (int ***)malloc(dims[0] * sizeof(int **));
    int ***data_recv_left = (int ***)malloc(dims[0] * sizeof(int **));
    int ***data_recv_right = (int ***)malloc(dims[0] * sizeof(int **));
    int ***data_recv_forward = (int ***)malloc(dims[0] * sizeof(int **));
    int ***data_recv_backward = (int ***)malloc(dims[0] * sizeof(int **));

    MPI_Request reqs[6];
    MPI_Status stats[6];

    MPI_Isend(&(data_send[0][0][0]), dims[0]*dims[1]*dims[2], MPI_INT, up_rank, 0, cart_comm, &reqs[0]);
    MPI_Irecv(&(data_recv_down[0][0][0]), dims[0]*dims[1]*dims[2], MPI_INT, down_rank, 0, cart_comm, &reqs[1]);
    MPI_Isend(&(data_send[0][0][0]), dims[0]*dims[1]*dims[2], MPI_INT, down_rank, 0, cart_comm, &reqs[2]);
    MPI_Irecv(&(data_recv_up[0][0][0]), dims[0]*dims[1]*dims[2], MPI_INT, up_rank, 0, cart_comm, &reqs[3]);
    MPI_Isend(&(data_send[0][0][0]), dims[0]*dims[1]*dims[2], MPI_INT, left_rank, 0, cart_comm, &reqs[4]);
    MPI_Irecv(&(data_recv_right[0][0][0]), dims[0]*dims[1]*dims[2], MPI_INT, right_rank, 0, cart_comm, &reqs[5]);
    MPI_Isend(&(data_send[0][0][0]), dims[0]*dims[1]*dims[2], MPI_INT, right_rank, 0, cart_comm, &reqs[6]);
    MPI_Irecv(&(data_recv_left[0][0][0]), dims[0]*dims[1]*dims[2], MPI_INT, left_rank, 0, cart_comm, &reqs[7]);
    MPI_Isend(&(data_send[0][0][0]), dims[0]*dims[1]*dims[2], MPI_INT, forward_rank, 0, cart_comm, &reqs[8]);
    MPI_Irecv(&(data_recv_backward[0][0][0]), dims[0]*dims[1]*dims[2], MPI_INT, backward_rank, 0, cart_comm, &reqs[9]);
    MPI_Isend(&(data_send[0][0][0]), dims[0]*dims[1]*dims[2], MPI_INT, backward_rank, 0, cart_comm, &reqs[10]);
    MPI_Irecv(&(data_recv_forward[0][0][0]), dims[0]*dims[1]*dims[2], MPI_INT, forward_rank, 0, cart_comm, &reqs[11]);

    MPI_Waitall(12, reqs, stats);

    // Print received data
    printf("Process %d: Received from up: %d, down: %d, left: %d, right: %d, forward: %d, backward: %d\n", rank,
           data_recv_up[0][0][0], data_recv_down[0][0][0], data_recv_left[0][0][0],
           data_recv_right[0][0][0], data_recv_forward[0][0][0], data_recv_backward[0][0][0]);

    // Free allocated memory
    for (int i = 0; i < dims[0]; ++i) {
        for (int j = 0; j < dims[1]; ++j) {
            free(data_send[i][j]);
        }
        free(data_send[i]);
    }
    free(data_send);

    for (int i = 0; i < dims[0]; ++i) {
        for (int j = 0; j < dims[1]; ++j) {
            free(data_recv_up[i][j]);
            free(data_recv_down[i][j]);
            free(data_recv_left[i][j]);
            free(data_recv_right[i][j]);
            free(data_recv_forward[i][j]);
            free(data_recv_backward[i][j]);
        }
        free(data_recv_up[i]);
        free(data_recv_down[i]);
        free(data_recv_left[i]);
        free(data_recv_right[i]);
        free(data_recv_forward[i]);
        free(data_recv_backward[i]);
    }
    free(data_recv_up);
    free(data_recv_down);
    free(data_recv_left);
    free(data_recv_right);
    free(data_recv_forward);
    free(data_recv_backward);

    MPI_Finalize();
    return 0;
}
``
