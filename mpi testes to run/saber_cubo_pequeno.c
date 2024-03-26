#include <stdio.h>

// Função para determinar em qual cubo menor um ponto está localizado
int determinarCubo(int x, int y, int z, int L, int lx, int ly, int lz) {
    int n_x = x / (L / lx);
    int n_y = y / (L / ly);
    int n_z = z / (L / lz);
    return n_x + n_y * (L / lx) + n_z * (L / ly) * (L / lz);
}

int main() {
    // Tamanho do cubo grande e do cubo menor
    int L, lx, ly, lz;

    // Obter o tamanho do cubo grande e do cubo menor
    printf("Digite o tamanho do cubo grande (L): ");
    scanf("%d", &L);
    
	printf("Digite o tamanho do cubo menor (l): ");
    scanf("%d", &lx);
	
    printf("Digite o tamanho do cubo menor (l): ");
    scanf("%d", &ly);

    printf("Digite o tamanho do cubo menor (l): ");
    scanf("%d", &ly);



    // Obter as coordenadas do ponto
    int x, y, z;
    printf("Digite as coordenadas do ponto (x y z): ");
    scanf("%d %d %d", &x, &y, &z);

    // Verificar em qual cubo menor o ponto está localizado
    int cubo = determinarCubo(x, y, z, L, lx, ly , lz);
    
    // Imprimir o resultado
    printf("O ponto (%d, %d, %d) está no cubo %d.\n", x, y, z, cubo);

    return 0;
}
