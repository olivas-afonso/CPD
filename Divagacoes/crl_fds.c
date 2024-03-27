#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>

int main (void){
    int n= 5;
    int p =4;
    int res;
    long long total=0;
    long long inicial=0;
    int y=0;

    int p_2= p*p;
    int n_2 = n*n;
    inicial = n*n;

    for (int x = 0; x < p; x++)
    {   
        res = floor (n_2/p_2);
        n_2 = n_2-res;
        p_2=p_2-1; 
        total = total + res; 
        printf ("x=%d\n",res);
        for (int z = 0; z < p && z<p_2; z++)
        {
            res = floor (n_2/p_2);
            n_2 = n_2-res;
            p_2=p_2-1; 
            total = total + res; 
            printf ("Res=%d\n",res);
        }     
    }
    
        
    printf ("Inicial:%lld Total:%lld\n",inicial, total);
        
}   

   /* int x_blocks = (NX + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int y_blocks = (NY + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int z_blocks = (NZ + BLOCK_SIZE - 1) / BLOCK_SIZE;
*/