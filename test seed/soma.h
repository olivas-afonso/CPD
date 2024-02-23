#ifndef SOMA_H
#define SOMA_H

#include <stdio.h>
#include <stdlib.h>


#define N_SPECIES 9

char *** grid;
unsigned int seed;
//long long cont_species[9]={0,0,0,0,0,0,0,0,0};
//long long cont1=0, cont2=0, cont3=0, cont4=0, cont5=0, cont6=0, cont7=0, cont8=0, cont9=0;
//long long max_cont1, max_cont2, max_cont3, max_cont4, max_cont5, max_cont6, max_cont7, max_cont8, max_cont9;
//long long max_count[9]={0,0,0,0,0,0,0,0,0};
//int max_gen[9];
//int current_gen =0;
//int max_gen1, max_gen2, max_gen3, max_gen4, max_gen5, max_gen6, max_gen7, max_gen8, max_gen9; 



void ver_celula_grid(int x , int y, int z);

void ver_grid(int N);


void ver_vizinhos();
int verif_regras_viva(int prox_estado);

int revive_celula(int valor_vizinhos[9]);
	
int verif_vizinhos( int pos_x, int pos_y ,int pos_z);

void geracao();

#endif
