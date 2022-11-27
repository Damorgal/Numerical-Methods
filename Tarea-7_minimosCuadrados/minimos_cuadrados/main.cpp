// Programa de minimos cuadrados
// Para compilar: c++ -o pro main.cpp funciones.cpp

// Tarea hecha por Diego Aaron Moreno Galvan

#include "cuadrados.hpp"

// La funcion cuadratica para el problema
double f(double* coef, double x)  {
    double ans;
    ans = coef[0]*x*x + coef[1]*x + coef[2];
    return ans;
}

int main(int narg, char **varg)  {
    int m = 50;

    // Creamos los valores de los puntos {xi, sin(xi)} para i = 1, ..., m
    double *x = crea_valores_xi(m);
    double *y = crea_valores_yi(x,m);

    double *alfa = minimos_cuadrados(x,y,m);

    // Guarda matriz en un archivo
    double **mat = crea_matriz(3,m);
    for(int i=0; i<m; i++)  {
        mat[0][i] = x[i];
        mat[1][i] = y[i];
        mat[2][i] = f(alfa,x[i]);
    }
    char s[] = "MinimosCuadrados.txt";
    guarda_matriz_txt(mat,3,m,s);

    // Imprimimos los coeficientes y el error
    printf("Los coeficientes de la parabola encontrados son:\n");
    printf("a = %lf,\nb = %lf,\nc = %lf\nParabola: %lfx^2 + %lfx + %lf\n",alfa[0],alfa[1],alfa[2],alfa[0],alfa[1],alfa[2]);
    double error = calcula_error_mincuadrados(mat[2],y,m);
    printf("Y el error E(xi) = 1/2SUMA( (f(xi) - yi)^2 ) = %lf\n",error);

    free(x);
    free(y);
    free(alfa);
    free(mat);
    return 0;
}
