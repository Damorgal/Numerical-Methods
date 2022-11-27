#include "interpolacion.hpp"

int main(int narg, char **varg)
{
    int nr,nc;
    // Leer los archivos de los puntos (xi ,yi)
    if(narg != 4) {cout << "Introduce la direccion del archivo del arreglo, le numero n y el lamda."<< endl; return 0;}
    double **A = lee_matriz(varg[1],&nr,&nc);
    double lamda = atof(varg[3]);
    int n = atoi(varg[2]);

    // Polinomio 
    double *puntosX = crea_vector(nr), *puntosY = crea_vector(nr);
    for(int i=0; i<nr; i++)	{
		puntosX[i] = A[i][0];
		puntosY[i] = A[i][1];
	}
	freeMat(A);
	int m = nr-1;
	
	// Discretizamos
	double a = 0, b = 10;
	double h = (b-a)/(double)n;
	double *x = crea_vector(n+1);
	for(int i=0; i<=n; i++)
		x[i] = a + h*(double)i;
		
	// Algoritmo
    double *fis = minimiza_diferencias(puntosX,puntosY,x,lamda,h,n,m);
	double **datos = evalua_datos_nubes(puntosX,fis,x,n,m);
	
	// Imprimimos la matriz de datos en txt para graficar
	char s[] = "puntos_nubes.txt";
	guarda_matriz_txt(datos,2,m+1,s);

    // Liberamos memoria
    freeMat(datos);
    free(fis);
    free(x);
    free(puntosX); free(puntosY);
    return 0;
}
