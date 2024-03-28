#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Definindo o tamanho do array
    const int N = 5;

    // Cada processo MPI tem seu próprio array de números inteiros
    int local_array[N];

    // Inicializando o array local com valores específicos para cada processo
    for (int i = 0; i < N; i++) {
        local_array[i] = (rank + 1) * (i + 1); // Por exemplo, processo 0: [1, 2, 3, 4, 5], processo 1: [2, 4, 6, 8, 10], etc.
		//printf("Array inicial: %d\n", local_array[i]);
	}

    // Array para armazenar o resultado da redução
    int global_array[N];

    // Reduzindo os arrays locais em um array global usando a operação de soma
    MPI_Reduce(local_array, global_array, N, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Somente o processo 0 imprime o array global resultante
    if (rank == 0) {
        printf("Array reduzido:\n");
        for (int i = 0; i < N; i++) {
            printf("%d ", global_array[i]);
        }
        printf("\n");
    }

    MPI_Finalize();
    return 0;
}
