#include "soma.h"

int vizinhos[26][3];
int vizinho_dezena = 0;
int valor;
int prox_estado = 0;
int valor_vizinhos[9];
		

void ver_celula_grid(int x , int y, int z){
		//("grid %c", grid[x][y][z]);	
}

/*
void ver_grid(int N){
	//("grid\n");
	for(int z = 0; z<N; z++){
		//("Layer %d\n", z);
		for(int y = 0; y<N; y++){
			for(int x = 0; x<N; x++){
				//("%d ", grid[x][y][z]);
			}
			//("\n");
		}
		//("\n");
	}
	//("-------------------\n");		
}
*/


void ver_grid(int N){
	//("grid\n");
	for(int x = 0; x<N; x++){
		//("Layer %d\n", x);
		for(int y = 0; y<N; y++){
			for(int z = 0; z<N; z++){
				//("%d ", grid[x][y][z]);
			}
			//("\n");
		}
		//("\n");
	}
	//("-------------------\n");	
	
	
}


void ver_vizinhos(){
	int cont = 0;
	int x = 0;
	int y = 0;
	int z = 0;
	int pos_val = 0;
	
		for( z = -1; z < 2; z++){
			for( y = -1; y < 2; y++){
				for( x = -1; x < 2; x++){
					if( y == 0) { 
						pos_val++;
					}
					if(z == 0) { 
						pos_val++;
					}
					if( x == 0){ 
						pos_val++;
					}
					cont++;		
					x = x;
					if(pos_val != 3) //("cont : %d pos = %d || z = %d || y = %d || x = %d \n", cont, pos_val, x ,y, z);
					pos_val = 0;
				}			
			}
		}
}
/*
void verifica_celula(){
	//("%d ", grid[x][y][z]);
	if( grid[x][y][z] 
}
*/

int verif_regras_viva(int prox_estado){
	int pos_x, pos_y, pos_z;
	int x, y, z;
	
	//("prox_estado:%d\n", prox_estado);
	
	if(prox_estado < 5){
		//celula morre
//		//("morre\n");
		return 0;
	}else{
		if(prox_estado < 14){
			//celula morre
//			//("vive\n");
			return 1;			
		}else{
			// celula morre
			//("morre\n");
			return 0;	
	
		}	
	}	
}

int revive_celula(int valor_vizinhos[9]){
	int i=0;
	int val_rep = 0;
	int n_rep = 0;
	
	for( i = 10; i > -1; i--){
		//("valor %d aparece %d\n", i, valor_vizinhos[i]);
		if(n_rep <= valor_vizinhos[i]){
			n_rep = valor_vizinhos[i];
			val_rep = i+1;
			//("maior val %d\n", val_rep);
		}
	}
	////("valor mais visto %d\n", maior_val);
return val_rep;
}

/*void apaga_array(){
	for(int contador = 0; contador < 11; contador++){
		valor_vizinhos[contador] = 0;	
	}
}
*/
void corre_vizinhos(int pos_x, int pos_y ,int pos_z){
	int x,y,z;
	int x_v,y_v,z_v;
	int pos_val, cont = 0;
	
	for( x = -1; x < 2; x++){
			for( y = -1; y < 2; y++){
				for( z = -1; z < 2; z++){
					
					/* salta a origem */
					if( y == 0) { 
						pos_val++;
						if(z == 0) { 
							pos_val++;
							if( x == 0){ 
								pos_val++;
							}
						}
					}	
						
					/* segue para o prox*/	
						
					if(pos_val != 3){
						
						x_v = pos_x+x;
						y_v = pos_y+y;
						z_v = pos_z+z;
						
						////("vizinho1: %d, %d, %d\n", x_v, y_v, z_v);
						
						if(x_v < 0){
							x_v = 127; 	//trocar o 4 por N
						}
						if(y_v < 0){
							y_v = 127;	//trocar o 4 por N
						}
						if(z_v < 0){
							z_v = 127;	//trocar o 4 por N
						}
						
						if(x_v > 3){
							x_v = 127; 	//trocar o 4 por N
						}
						if(y_v > 3){
							y_v = 127;	//trocar o 4 por N
						}
						if(z_v > 3){
							z_v = 127;	//trocar o 4 por N
						}
						
						////("vizinho: %d, %d, %d\n",x_v, y_v, z_v);
						
						vizinhos[cont][0] = x_v;
						vizinhos[cont][1] = y_v;
						vizinhos[cont][2] = z_v;
						cont++;	
						
						if(grid[x_v][y_v][z_v] > 9){
							vizinho_dezena = (grid[x_v][y_v][z_v] / 10); 
							}else{
							vizinho_dezena = grid[x_v][y_v][z_v];
						}
						////("valor na grid %d e outro %d\n",  grid[x_v][y_v][z_v], vizinho_dezena );
						
						if(vizinho_dezena > 0){
							prox_estado++;
							valor = vizinho_dezena - 1;
							//("vetor valor:%d\n", valor);
							valor_vizinhos[valor]++;
						}
						
						////("prox _ estado:%d\n", prox_estado);
						
						//prox_estado += verifica_vivo(x, y, z) + verifica_n_vivo(x, y, z);
					
					}
					
					pos_val = 0;
				}			
			}
		}
}


int verif_vizinhos(int pos_x, int pos_y ,int pos_z){
	
	int ret = 0;
	
	//apaga_array(valor_vizinhos);
	for(int contador = 0; contador < 11; contador++){
		valor_vizinhos[contador] = 0;	
	}
	
	corre_vizinhos(pos_x, pos_y, pos_z);
			
	if(grid[pos_x][pos_y][pos_z] == 0){
		////("celula_morta: %d\n", prox_estado);
		if(prox_estado < 11){
			if(prox_estado >6){
//				//("vive\n");
				grid[pos_x][pos_y][pos_z] = revive_celula(valor_vizinhos) * 10;	
			}else{
				//grid[pos_x][pos_y][pos_z] = 0 ;
				return 0;
			}
		}
	}else{
		////("valor:%d\n", grid[pos_x][pos_y][pos_z]);
		////("valor prev %d\n", vizinho_dezena);
		ret = verif_regras_viva(prox_estado); 
		if( ret != 1){
			//deixa a grid nas unidades
			//grid[pos_x][pos_y][pos_z] = 0;
		}else{
		//	//("valor desta gen :%d\n", grid[pos_x][pos_y][pos_z] );
			grid[pos_x][pos_y][pos_z] += grid[pos_x][pos_y][pos_z] *10;	
		//	//("valor prox gen :%d\n", grid[pos_x][pos_y][pos_z]);
			
		}

	}


}

void geracao(){
	int cont=0;
	for( int x = 0; x < 128; x++){
		//("\nproximo layer: %d\n\n\n", x);
		
		for( int y = 0; y < 128; y++){
			for( int z = 0; z < 128; z++){
				cont++;
				//("\nproxima celula: %d\nvalor: %d\n", cont, grid[x][y][z]);
				verif_vizinhos(x,y,z);
				prox_estado = 0;
			}
		}
	}	
}
