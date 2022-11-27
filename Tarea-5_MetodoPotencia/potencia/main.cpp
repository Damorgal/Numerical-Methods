#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

double* potencia(double** A, int omax, double tol, int n);
double producto_punto(double* a, double* b, int n);
void producto_MatVec(double** P, double* a, double* b, int n);
double error_eigenpar(double *y, double lm, double *v, int n);
double norma(double *vec, int n);

void inicializa_vec(double* vec, double val, int n);
double **lee_matriz(char *dir, int *nr, int *nc);
double* crea_vector(int n);
double** crea_matriz(int nr, int nc);
void freeMat(double** mat);

int main(int narg, char **varg)  {
    int nr,nc,nmax;
    double tolerancia;

    if(narg != 4) {printf("\nIntroduce la direccion del archivo binario, numero maximo de operaciones y tolerancia\n"); return 0;}
    char *dirmat = varg[1];
    nmax = atoi(varg[2]); tolerancia = atof(varg[3]);

    // Leemos los archivos
    double **A = lee_matriz(dirmat,&nr,&nc);
    if(nr != nc) {printf("Error, las matrices en el programa deben ser cuadradas.\n"); freeMat(A); return 0;}
    int n = nr;

    // Hacemos el metodo de la potencia para encontrar el eigenpar tal que el va.p. |lamda| es el mayor
    // La primera entrada de eigpar es el va.p. mas grande lamda y las siguientes corresponden al eigenvector
    double *eigpar = potencia(A, nmax, tolerancia, n);

    free(eigpar);
    freeMat(A);
	return 0;
}

double* potencia(double** A, int omax, double tol, int n)   {
    // Inicializamos vector v0;
	double *v = crea_vector(n);
	inicializa_vec(v,1,n);

	// Creamos vector y y lo inicializamos
	double *y = crea_vector(n);
	int k=0; double lm=0;
	producto_MatVec(A,v,y,n);

	// Hacemos algoritmo de la potencia
	double *aux = crea_vector(n);
	for(; k < omax && error_eigenpar(y,lm,v,n) >= tol; k++)  {
        // Hacemos operacion y = Av
        if(k > 0)
            producto_MatVec(A,v,y,n);

        // Hacemos operacion v = y/||y||
        double nor = norma(y,n);
        for(int i=0; i<n; i++)
            v[i] = y[i] / nor;

        // Hacemos operacion lamda = vAv
        producto_MatVec(A,v,aux,n);
        lm = producto_punto(v,aux,n);
	}

	printf("El tamagno de la matriz es de: %dx%d\n",n,n);
	printf("El valor propio mas grande encontrado es: %lf\n",lm);
	printf("El numero maximo de iteraciones era de: %d y las realizadas fueron: %d\n",omax,k);
	printf("El error calculado de |Av-lv| es: %lf (tolerancia = %lf)\n",error_eigenpar(y,lm,v,n),tol);
    double *x = crea_vector(n+1);
    x[0] = lm; for(int i=0; i<n; i++) x[i+1] = v[i];
    free(v); free(y); free(aux);
    return x;
}

double norma(double *vec, int n)   {
    double sum=0;
    for(int i=0; i<n; i++)
        sum += vec[i]*vec[i];
    sum = sqrt(sum);
    return sum;
}

// Funcion que hace el producto punto de los vectores a y b (le debes pasar la longitud n de los vectores)
double producto_punto(double* a, double* b, int n)    {
    double pp=0;
	#pragma omp parallel for reduction(+:pp)
	for (int i=0; i<n; i++)
        pp += a[i] * b[i];

    return pp;
}

// Funcion que hace el producto de una matriz con un vector (P y a) y lo guarda en b (le debes pasar la longitud n de los vectores)
void producto_MatVec(double** P, double* a, double* b, int n)   {
    #pragma omp parallel for
	for (int i=0; i<n; i++) {
        double sum=0;
        //#pragma omp parallel for reduction(+:sum)
        for(int j=0; j<n; j++)
            sum += P[i][j] * a[j];
        b[i] = sum;
	}
}

// Calcula el error de que una matriz por un vector x de la solucion del vector b
double error_eigenpar(double *y, double lm, double *v, int n)  {
    double *aux = crea_vector(n);
    for(int i=0; i<n; i++)
        aux[i] = y[i]-(lm * v[i]);
    double nor = norma(aux,n);
    free(aux);
    return nor;
}


// Inicializa todos las entradas de una matriz a un valor dado
void inicializa_vec(double* vec, double val, int n)   {
    for(int i=0; i<n; i++)
            vec[i] = val;
}

// Funcion que lee de archivo una matriz
double **lee_matriz(char *dir, int *nr, int *nc) {
    double **mat;

    FILE *fp = fopen(dir, "rb");
    if (fp==NULL) {printf("El archivo %s no puedo ser leido.\n", dir); exit (1);}

    // Lectura del tama~no de la matriz
    fread(nr, sizeof(int), 1, fp);
    fread(nc, sizeof(int), 1, fp);
    // Reservamos memoria
    mat = crea_matriz(*nr, *nc);
    // Lectura de los datos
    fread(mat[0], sizeof(double), (*nr)*(*nc), fp);
    fclose(fp);
    return(mat);
}

double* crea_vector(int n)  {
    double* aux;
    aux = (double*)malloc(sizeof(double)*n);
    return aux;
}

double** crea_matriz(int nr, int nc) {
    int i;
    double** mat;

    mat = (double**)malloc(nr * sizeof(double*));
    if(mat==NULL) return(NULL);
    mat[0] = (double*)malloc(nr * nc * sizeof(double));
    if(mat[0]==NULL) return(NULL);

    for(i=1; i<nr; i++)
        mat[i] = mat[i-1] + nc;

    return mat;
}

//libera memoria de una matriz
void freeMat(double** mat)   {
    free(mat[0]);
    free(mat);
}
