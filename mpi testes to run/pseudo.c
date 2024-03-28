sub_div_x[2] = [2, 2, 1];
sub_div_y [2] = [2, 3];
sub_div_z [2] = [2, 3];

int limite_inf_x;
int limite_inf_y;
int limite_inf_z;


void limites_x (){
    if (coodernada de x do pc == 0){
        limite_inf_x = 0; 
    }
        
    else{
        for (int i = 0; i < coodernada de x do pc -1; ++ i){
            limite_inf_x = limite_inf_x + sub_div_x [i];
        }
    }

    for (int i = 0; i< coodernada de y do pc; ++i){
        limite_sup_x = limite_sup_x + sub_div_x [i];
    }
}

void limites_y (){
    if (coodernada de y do pc == 0){
        limite_inf_y = 0; 
    }
        
    else{
        for (int i = 0; i < coodernada de y do pc -1; ++ i){
            limite_inf_y = limite_inf_y + sub_div_y [i];
        }
    }

    for (int i = 0; i< coodernada de y do pc; ++i){
        limite_sup_y = limite_sup_y + sub_div_y [i];
    }
}

void limites_z (){
    if (coodernada de z do pc == 0)
        limite_inf_z = 0;
    else{
        for (int i = 0; i < coodernada de z do pc -1; ++ i){
            limite_inf_z = limite_inf_z + sub_div_z [i];
        }
    }

    for (int i = 0; i< coodernada de z do pc; ++i){
        limite_sup_z = limite_sup_z + sub_div_z [i];
    }
}


int main

int varrimento_x = 1;
int varrimento_y = 1;
int varrimento_z = 1;

limites_x ();
limites_y ();
limites_z();

for (x...){
    if (x>=limite_inf_x && x<limite_sup_x){
        flag_y = 1;
        ++varrimento_x;
    }
    
    for (y...){
        if (y>=limite_inf_y && y<limite_sup_y){
            flag_y = 1;
            ++varrimento_y;
        }

        for (z...){
            if (z>=limite_inf_y && z<limite_sup_y && flag_x = 1 && flag_y == 1){
                 matriz [varrimento_x-1][varrimento_y-1][varrimento_z] ;
                 ++varrimento_z;
            }
        }
        varrimento_z = 0;
    }
    varrimento_y = 0;
}
    