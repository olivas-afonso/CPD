#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Define Cartesian Topology
    int dims[2] = {2, 2};  // 2x2 grid
    int periods[2] = {1, 1};  // Enable wraparound
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);

    // Get Cartesian coordinates of current process
    int my_coords[2];
    MPI_Cart_coords(cart_comm, rank, 2, my_coords);

    // Get neighbors
    int up_rank, down_rank, left_rank, right_rank;
    MPI_Cart_shift(cart_comm, 0, 1, &up_rank, &down_rank);
    MPI_Cart_shift(cart_comm, 1, 1, &left_rank, &right_rank);

    // Send and receive data between neighbors
    int data_send = rank * 10 + my_coords[0] * 100 + my_coords[1];  // Example data
    int data_recv_up, data_recv_down, data_recv_left, data_recv_right;
    MPI_Sendrecv(&data_send, 1, MPI_INT, up_rank, 0, &data_recv_down, 1, MPI_INT, down_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send, 1, MPI_INT, down_rank, 0, &data_recv_up, 1, MPI_INT, up_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send, 1, MPI_INT, left_rank, 0, &data_recv_right, 1, MPI_INT, right_rank, 0, cart_comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&data_send, 1, MPI_INT, right_rank, 0, &data_recv_left, 1, MPI_INT, left_rank, 0, cart_comm, MPI_STATUS_IGNORE);

    // Print received data
    printf("Process %d: Received from up: %d, down: %d, left: %d, right: %d\n", rank, data_recv_up, data_recv_down, data_recv_left, data_recv_right);

    MPI_Finalize();
    return 0;
}
