#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int factoriza_LU(double** A, double** L, double** U, int n);
double* resuelve_triangularL(double** A, double* b, int n, int *ind);
double* resuelve_triangularU(double** A, double* b, int n, int *ind);
void calcula_errorM(double **A, double *b, double *x, int n);
void calcula_errorLU(double **A, double **L, double **U, int n);

void inicializa(double** mat, double val, int nr, int nc);
double *lee_arreglo(char *dir, int *n);
double **lee_matriz(char *dir, int *nr, int *nc);
double* crea_vector(int n);
double** crea_matriz(int nr, int nc);
void freeMat(double** mat);

int main(int narg, char **varg)  {
    int n,nr,nc;
    double *y=NULL,*x=NULL;
    char *dirmat = varg[1];
    char *dirvec = varg[2];

    if(narg != 3) {printf("\nIntroduce la direccion de los DOS archivos binarios (matriz y sol. b)\n"); return 0;}
    // Leemos los archivos
    double **A = lee_matriz(dirmat,&nr,&nc);
    if(nr != nc) {printf("Error, las matrices en el programa deben ser cuadradas.\n"); freeMat(A); return 0;}
    double *b = lee_arreglo(dirvec, &n);
    if(n != nc) {printf("Error, la dimension del vector b y la matriz no coinciden.\n"); freeMat(A); free(b); return 0;}

	// Encontramos factorizacion LU
	double **L = crea_matriz(n,n);
    double **U = crea_matriz(n,n);
    // Indicador ind nos dira si se pudo factorizar A como LU
    int ind = factoriza_LU(A,L,U,n);

    // Resolvemos el sistema
    // indicador inf para ver si tiene infinitas soluciones
    int inf=0;
    if(ind == 0)    {
        y = resuelve_triangularL(L,b,n,&inf);
        if(y != NULL && inf == 0)
            x = resuelve_triangularU(U,y,n,&inf);
    }

    // Salida del programa
    if(ind == 1) printf("No se puede realizar la factorizacion LU.\n");
    else {
        printf("La matriz SI tiene factorizacion LU\n");
        calcula_errorLU(A,L,U,n);

        if(inf == 1) printf("La matriz tiene infinitas soluciones.\n");
        else if(x == NULL) printf("La matriz no tiene solucion.\n");
        else    {
            printf("\nEl tamagno de la matriz es: %dx%d\n", nr, nc);
            printf("La matriz SI tiene solucion\n");
            calcula_errorM(A,b,x,n);
        }
    }
    free(y);
	free(x);
    free(b);
    freeMat(L);
    freeMat(U);
    freeMat(A);
	return 0;
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

// Calcula el error de que una matriz por un vector x de la solucion del vector b
void calcula_errorM(double **A, double *b, double *x, int n)  {
    double sum=0;
    for(int i=0; i<n; i++)  {
        double aux=0;
        for(int j=0; j<n; j++)
            aux += A[i][j]*x[j];
        sum += std::pow(aux - b[i], 2);
    }
    sum = std::sqrt(sum);
    printf("Y el error calculado de la solucion es: %lf\n", sum);
}

// Calcula el error de que una matriz A sea igual a dos matrices L y U multiplicadas (Norma de Frobenius)
void calcula_errorLU(double **A, double **L, double **U, int n) {
    double sum=0;
    for(int i=0; i<n; i++)  {
        for(int j=0; j<n; j++)  {
            double aux=0;
            for(int k=0; k<n; k++)
                aux += L[i][k]*U[k][j];
            sum += std::pow(A[i][j] - aux, 2);
        }
    }
    sum = std::sqrt(sum);
    printf("Y el error calculado de la matriz LU es: %lf\n", sum);
}

// Inicializa todos las entradas de una matriz a un valor dado
void inicializa(double** mat, double val, int nr, int nc)   {
    for(int i=0; i<nr; i++)
        for(int j=0; j<nc; j++)
            mat[i][j] = val;
}

// Funcion que lee de archivo un vector
double *lee_arreglo(char *dir, int *n) {
    FILE *fp;
    double *vec;

	fp = fopen(dir, "rb");
	if (fp==NULL) {printf("El archivo %s no puedo ser leido.\n", dir); exit (1);}

	// Pedir memoria
	fread(n, sizeof(int), 1, fp);
    vec = crea_vector(*n);

	// Leer arreglo
	fread(vec, sizeof(double), *n, fp);

	fclose(fp);
	return(vec);
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
