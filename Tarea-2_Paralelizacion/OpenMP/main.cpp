#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

void suma(double* a, double* b, double* c, int n);
void multiplica(double* a, double* b, double* c, int n);
double producto_punto(double* a, double* b, int n);
void producto_MatVec(double** P, double* a, double* b, int n);
void producto_matricial(double** P, double** Q, double** R, int n);

double* crea_vector(int n);
double** crea_matriz(int nr, int nc);
void inicializa_vector(double *vec, int n);
void inicializa_matriz(double **A, int n);
void freeMat(double** mat);
void tiempoTranscurrido(struct timeval t1, struct timeval  t2);

int main(int narg, char *varg[])  {
    double ppunto, *a, *b, *c, **P, **Q, **R;
    struct timeval t1, t2;
    int n, threads;

    int m = omp_get_max_threads();
	printf("omp_get_max_threads = %i\n", m);

    // Debes correr el programa por ejemplo como: ./programa 100000 4
    if(narg != 3) {printf("\nIntroduce el tamagno del vector y el numero de threads\n"); exit(1);}
    n = atoi(varg[1]);
    threads = atoi(varg[2]);

    a = crea_vector(n); inicializa_vector(a,n);
	b = crea_vector(n); inicializa_vector(b,n);
	c = crea_vector(n);
	P = crea_matriz(n,n); inicializa_matriz(P,n);
	Q = crea_matriz(n,n); inicializa_matriz(Q,n);
    R = crea_matriz(n,n);

    // Iniciamos con el numero de threads dados
    omp_set_num_threads(threads);

    // Tomamos el tiempo que le cuesta hacer la suma
    gettimeofday(&t1, NULL);
    suma(a, b, c, n);
    gettimeofday(&t2, NULL);
    printf("\nEl tiempo que toma hacer la suma:\n");
    tiempoTranscurrido(t1, t2);

    // Tomamos el tiempo que le cuesta hacer la multiplicacion
    gettimeofday(&t1, NULL);
    multiplica(a, b, c, n);
    gettimeofday(&t2, NULL);
    printf("\nEl tiempo que toma hacer la multiplicacion:\n");
    tiempoTranscurrido(t1, t2);

	free(c);
	// Tomamos el tiempo que le cuesta hacer el producto punto
    gettimeofday(&t1, NULL);
    ppunto = producto_punto(a, b, n);
    gettimeofday(&t2, NULL);
    printf("\nEl producto punto es: %lf\nEl tiempo que toma hacer el producto punto:\n",ppunto);
    tiempoTranscurrido(t1, t2);

    // Tomamos el tiempo que le cuesta hacer el producto de una matriz con un vector
    gettimeofday(&t1, NULL);
    producto_MatVec(P, a, b, n);
    gettimeofday(&t2, NULL);
    printf("\nEl tiempo que toma hacer el producto de una matriz con un vector:\n");
    tiempoTranscurrido(t1, t2);

    free(b);
    free(a);
    // Tomamos el tiempo que le cuesta hacer el producto de dos matrices
    gettimeofday(&t1, NULL);
    producto_matricial(P, Q, R, n);
    gettimeofday(&t2, NULL);
    printf("\nEl tiempo que toma hacer el producto de dos matrices:\n");
    tiempoTranscurrido(t1, t2);

    freeMat(P);
    freeMat(Q);
    freeMat(R);
	return 0;
}

// Funcion que suma los vectores a y b y los guarda en c (le debes pasar la longitud n de los vectores)
void suma(double *a, double *b, double *c, int n) {
	#pragma omp parallel for
	for (int i=0; i<n; i++)
        c[i] = a[i] + b[i];
}

// Funcion que multiplica los vectores a y b y los guarda en c (le debes pasar la longitud n de los vectores)
void multiplica(double* a, double* b, double* c, int n)   {
    #pragma omp parallel for
	for (int i=0; i<n; i++)
        c[i] = a[i] * b[i];
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

// Funcion que hace el producto de dos matrices P y Q y lo guarda en R (le debes pasar la longitud n de los vectores)
void producto_matricial(double** P, double** Q, double** R, int n)  {
    #pragma omp parallel for
	for (int i=0; i<n; i++) {
        //#pragma omp parallel for
        for(int j=0; j<n; j++)  {
            double sum=0;
            //#pragma omp parallel for reduction(+:sum)
            for(int k=0; k<n; k++)
                sum += P[i][k] * Q[k][j];
            R[i][j] = sum;
        }
	}
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

void inicializa_vector(double *vec, int n)  {
    for(int i=0; i<n; i++)
        vec[i] = ((rand()/(double)RAND_MAX)*2)-1;
}

void inicializa_matriz(double **A, int n)  {
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            A[i][j] = ((rand()/(double)RAND_MAX)*2)-1;
}

void freeMat(double** mat)   {
    free(mat[0]);
    free(mat);
}

void tiempoTranscurrido(struct timeval t1, struct timeval  t2) {
    double elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;  // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;      // us to ms
    printf("Tiempo transcurrido %6.4f mseg. \n", elapsedTime);
}
