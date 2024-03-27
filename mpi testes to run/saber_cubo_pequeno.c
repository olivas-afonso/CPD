#include <stdio.h>

// Função para determinar em qual cubo menor um ponto está localizado
int determinarCubo(int x, int y, int z, int L, int lx, int ly, int lz) {
    int n_x = x / (L / lx);
    int n_y = y / (L / ly);
    int n_z = z / (L / lz);
    return n_x + n_y * (L / lx) + n_z * (L / ly) * (L / lz);
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

int main() {
    // Tamanho do cubo grande e do cubo menor
    int L, l, lx, ly, lz;

    // Obter o tamanho do cubo grande e do cubo menor
    printf("Digite o tamanho do cubo grande (L): ");
    scanf("%d", &L);
    
	printf("Digite o tamanho do cubo menor (l): ");
    scanf("%d", &l);
	
	
	 int maior_linha = 0, maior_linha_prev = L ;
    int a, b, c;
	int a_final, b_final, c_final;
	

    for (a = 1; a <= l; a++) {
        if (x % a == 0) {
            for (b = a; b <= l / a; b++) {
                if (l % (a * b) == 0) {
                    c = l / (a * b);
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

	
	
	
	// Verificar se y é maior que 0
    if (l <= 0) {
        printf("Erro: O tamanho do array deve ser maior que 0.\n");
        return 1; // Terminar o programa com erro
    }

    // Criar um array com o tamanho especificado
    int array[l];
	
	divide_number_parts(L, l, array);
	
    // Imprimir o array
    printf("Array gerado:\n");
    for (int i = 0; i < y; i++) {
        printf("%d ", array[i]);
    }
    printf("\n");

	
	lx = array[0];	
	ly = array[1];		
	lz = array[2];	


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
