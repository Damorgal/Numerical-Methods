#include "interpolacion.hpp"

int main(int narg, char **varg)
{
    int nr,nc;
    // Leer los archivos de los puntos (xi ,yi)
    if(narg != 3) {cout << "Introduce la direccion del archivo del arreglo, y el numero m."<< endl; return 0;}
    double **A = lee_matriz(varg[1],&nr,&nc);
    int m = atoi(varg[2]);

    // Arreglos
    double *puntosX = crea_vector(nr), *puntosY = crea_vector(nr);
    for(int i=0; i<nr; i++)	{
		puntosX[i] = A[i][0];
		puntosY[i] = A[i][1];
	}
	freeMat(A);
	int n = nr-1;
		
	// Algoritmo
    double *M = coef_spline_cubico(puntosX,puntosY,n);
	double **datos = evalua_datos_spline(puntosX,puntosY,M,n,m);
	
	// Imprimimos la matriz de datos en txt para graficar
	char s[] = "puntos_splines2.txt";
	guarda_matriz_txt(datos,2,m+1,s);

    // Liberamos memoria
    freeMat(datos);
    free(M);
    free(puntosX); free(puntosY);
    return 0;
}
