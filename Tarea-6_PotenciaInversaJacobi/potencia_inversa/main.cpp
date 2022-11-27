#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_math.h>


typedef struct datos_potencia_inversa {
    double mu; // Valor propio
    double *x; // Vector propio
    int k;     // Numero de operaciones realizadas
    double e;  // Norma ||w - px^|| < tolerancia

} eigpar_pot_inv;

void valores_propios(double **A, double d, int N, int M, double tol, int n);
eigpar_pot_inv potencia_inversa(double** A, double *x, double del, int omax, double tol, int n);
double* resuelve_ecuacion(double **A, double *b, int n);
double producto_punto(double* a, double* b, int n);
double** resta_matricial_diagonal(double** A, double del, int n);
int factoriza_LU(double** A, double** L, double** U, int n);
double* resuelve_triangularL(double** A, double* b, int n, int *ind);
double* resuelve_triangularU(double** A, double* b, int n, int *ind);
double norma(double *vec, int n);
double norma_inf(double **A, int n);

void inicializa(double** mat, double val, int nr, int nc);
void inicializa_vec(double* vec, double val, int n);
double **lee_matriz(char *dir, int *nr, int *nc);
double* crea_vector(int n);
double** crea_matriz(int nr, int nc);
void freeMat(double** mat);

int main(int narg, char **varg)  {
    int nr,nc;

    if(narg != 3) {printf("\nIntroduce la direccion del archivo binario y el entero N.\n"); return 0;}
    char *dirmat = varg[1];
    int N = atoi(varg[2]);

    // Leemos los archivos
    double **A = lee_matriz(dirmat,&nr,&nc);
    if(nr != nc) {printf("Error, las matrices en el programa deben ser cuadradas.\n"); freeMat(A); return 0;}
    int n = nr;

    // Calculamos d = norma infinito de matriz
    double d = norma_inf(A,n);

	/****** Definiciones ******/
    int M = 1000;    // Numero maximo de iteraciones
    double rem = sqrt(DBL_EPSILON); // tolerancia de raiz de epsilon maquina
    double tolerancia = rem;

    // Funcion para encontrar los valores propios con una discretizacion de N intervalos
	valores_propios(A, d, N, M, tolerancia, n);

    freeMat(A);
	return 0;
}

void valores_propios(double **A, double d, int N, int M, double tol, int n)    {
    eigpar_pot_inv data;
	double *x = crea_vector(n);
    double mu0 = -10*d;
    double dta; // Se convergera al va.p. mas cercano a delta
    for(int t=0; t<=N; t++) {
        // Inicializamos vector v0;
        inicializa_vec(x,1,n);

        // Inicializamos el delta de la discretizacion
        dta = (2*d/N)*t - d;

        // Llamamos la funcion potencia inversa
        data = potencia_inversa(A, x, dta, M, tol, n);

        // Verificamos que e<tol y |mu - mu0|> 0.0001 e imprimimos
        if(data.e < tol && abs(data.mu - mu0) > 0.0001) {
            printf("Un valor propio encontrado es: %lf\n", data.mu);
            printf("El numero maximo de iteraciones era de: %d y las realizadas fueron: %d\n", M, data.k);
            printf("El error calculado de e=||w-px^|| es: %lf (tolerancia = %lf)\n", data.e, tol);
            mu0 = data.mu;
        }
    }
    free(x);
}

eigpar_pot_inv potencia_inversa(double** A, double *x, double del, int omax, double tol, int n)   {
	// Creamos vectores y variables
	double *y = crea_vector(n);
	double *aux = crea_vector(n);
	double *w = crea_vector(n);
	double *r = crea_vector(n);
	int k=0; double mu; double p;

	// Hacemos algoritmo de la potencia potencia inversa
	do  {
        // Resolvemos (A - del*I)y = x
        double **resta = resta_matricial_diagonal(A,del,n);
        y = resuelve_ecuacion(resta, x, n);
        freeMat(resta);

        // Hacemos operacion aux = y/||y||, y la operacion w = x/||y||
        double nor = norma(y,n);
        for(int i=0; i<n; i++)  {
            aux[i] = y[i] / nor;
            w[i] = x[i] / nor;
        }

        // Hacemos operacion p = aux * w
        p = producto_punto(aux,w,n);

        // Hacemos operacion mu = del + p
        mu = del + p;

        // Hacemos la operacion r = w - p*aux
        for(int i=0; i<n; i++)
            r[i] = w[i] - (p*aux[i]);

        // Hacemos la operacion x = aux
        for(int i=0; i<n; i++)
            x[i] = aux[i];

        k++;
	}
	while( k < omax && norma(r,n) >= tol);

    eigpar_pot_inv data;
    data.e = norma(r,n);
    data.k = k;
    data.mu = mu;
    data.x = x;

    free(w); free(r); free(y); free(aux);
    return data;
}

double* resuelve_ecuacion(double **A, double *b, int n) {
    // Encontramos factorizacion LU
	double **L = crea_matriz(n,n);
    double **U = crea_matriz(n,n);
    // Indicador ind nos dira si se pudo factorizar A como LU
    int ind = factoriza_LU(A,L,U,n);

    // Resolvemos el sistema
    // indicador inf para ver si tiene infinitas soluciones
    int inf=0;
    double *y = crea_vector(n);
    double *x = crea_vector(n);
    if(ind == 0)    {
        y = resuelve_triangularL(L,b,n,&inf);
        if(y != NULL && inf == 0)
            x = resuelve_triangularU(U,y,n,&inf);
    }

    // Salida del programa
    if(ind == 1) printf("No se puede realizar la factorizacion LU.\n");
    else {
        if(inf == 1) printf("La matriz tiene infinitas soluciones.\n");
        else if(x == NULL) printf("La matriz no tiene solucion.\n");
    }
    free(y);
    freeMat(L);
    freeMat(U);
    return x;
}

// Calcula la norma 2 de un vector
double norma(double *vec, int n)   {
    double sum=0;
    for(int i=0; i<n; i++)
        sum += vec[i]*vec[i];
    sum = sqrt(sum);
    return sum;
}

// Calcula la norma infinito de una matriz = max (suma de cada fila)
double norma_inf(double **A, int n) {
    double sum, mx;
    for(int i=0; i<n; i++)  {
        sum = 0;
        for(int j=0; j<n; j++)
            sum += abs(A[i][j]);
        if(i == 0) mx = sum;
        else if(mx < sum) mx = sum;
    }
    return mx;
}

// Factorizacion por Doolittle, regresa 0 si pudo factorizar, sino 1.
int factoriza_LU(double** A, double** L, double** U, int n)   {
    if(A[0][0] == 0) return 1;
    // Llenamos las matrices L y U de ceros
    inicializa(L,0,n,n);
    inicializa(U,0,n,n);

    for(int i=0; i<n; i++) {
        U[0][i] = A[0][i];
        L[i][0] = A[i][0] / U[0][0];
    }
    for(int i=1; i<n; i++)
        for(int j=i; j<n; j++) {
            // Calculamos el renglon i de U
            double sum=0;
            for(int k=0; k<i; k++)
                sum += L[i][k] * U[k][j];
            U[i][j] = A[i][j] - sum;
            // Calculamos la columna i de L
            if(U[i][i] == 0) return 1;
            sum=0;
            for(int k=0; k<i; k++)
                sum += L[j][k] * U[k][i];
            L[j][i] = (A[j][i] - sum) / U[i][i];
        }
    return 0;
}

// Resuelve una matriz triangular inferior
double* resuelve_triangularL(double** A, double* b, int n, int *ind)   {
    double *x = crea_vector(n);
    for(int i=0; i<n; i++)  {
        double sum=0;
        if(i == 0)  {
            if(A[i][i] != 0)
                x[i] = b[i]/A[i][i];
            else if(b[i] == 0)
                *ind=1;
            else {free(x); return NULL;}
        }
        else    {
            for(int j=0; j<i; j++)
                sum += A[i][j]*x[j];
            if(A[i][i] != 0)
                x[i] = (b[i] - sum)/A[i][i];
            else if((b[i]-sum) == 0)
                *ind=1;
            else {free(x); return NULL;}
        }
    }
    return x;
}

// Resuelve una matriz triangular superior
double* resuelve_triangularU(double** A, double* b, int n, int *ind)   {
    double *x = crea_vector(n);
    for(int i=n-1; i>=0; i--)  {
        double sum=0;
        if(i == n-1)  {
            if(A[i][i] != 0)
                x[i] = b[i]/A[i][i];
            else if(b[i] == 0)
                *ind=1;
            else {free(x); return NULL;}
        }
        else    {
            for(int j=i+1; j<n; j++)
                sum += A[i][j]*x[j];
            if(A[i][i] != 0)
                x[i] = (b[i] - sum)/A[i][i];
            else if((b[i]-sum) == 0)
                *ind=1;
            else {free(x); return NULL;}
        }
    }
    return x;
}

// Funcion que hace el producto punto de los vectores a y b (le debes pasar la longitud n de los vectores)
double producto_punto(double* a, double* b, int n)    {
    double pp=0;
	#pragma omp parallel for reduction(+:pp)
	for (int i=0; i<n; i++)
        pp += a[i] * b[i];

    return pp;
}

// Calcula la resta de A - del*I
double** resta_matricial_diagonal(double** A, double del, int n)    {
    double **aux = crea_matriz(n,n);
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)  {
            aux[i][j] = A[i][j];
            if(i == j)
                aux[i][j] -= del;
        }
    return aux;
}

// Inicializa todos las entradas de una matriz a un valor dado
void inicializa(double** mat, double val, int nr, int nc)   {
    for(int i=0; i<nr; i++)
        for(int j=0; j<nc; j++)
            mat[i][j] = val;
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
