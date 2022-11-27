#include "interpolacion.hpp"

int main(int narg, char **varg)
{
    int nr,nc,n;
    // Leer los archivos de la matriz A y el vector b
    if(narg != 2) {cout << "Introduce la direccion del archivo del arreglo."<< endl; return 0;}
    double **A= lee_matriz(varg[1],&nr,&nc);

    // Polinomio interpolador
    double *abscisas = crea_vector(nr), *ordenadas = crea_vector(nr);
    for(int i=0; i<nr; i++)	{
		abscisas[i] = A[i][0];
		ordenadas[i] = A[i][1];
	}
	freeMat(A);
	n = nr-1;
    double *c = coef_pol_interpolador(abscisas,ordenadas,n);
    double **datos = evalua_datos(abscisas,c,n);

	// Imprimimos la matriz de datos en txt para graficar
	char s[] = "puntos_interpolacion_n14.txt";
	guarda_matriz_txt(datos,2,4*n+1,s);
	
    // Liberamos memoria
    freeMat(datos);
    free(c);
    free(abscisas);
    free(ordenadas);
    return 0;
}
