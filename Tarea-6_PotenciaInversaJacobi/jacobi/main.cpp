#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_math.h>


typedef struct datos_jacobi {
    double *lm; // Valores propios de matriz
    double **V; // Vectores propios de matriz
    int k;     // Numero de operaciones realizadas
    double e;  // Entrada mas grande de matriz A^(k-1)

} jacobiData;

typedef struct indices_matriz {
    int i;
    int j;
    double mx; // Los indices del valor maximo

} indices;

jacobiData jacobi(double **A, int omax, double tol, int n);
void calcula_error_jacobi(double **A, double **V, double *lm, int n);
double** traspuesta(double **A, int n);
double** rotacion_givens(double cos, double sen, int coori, int coorj, int n);
indices busca_mayor(double **A, int n);

void producto_matricial(double** P, double** Q, double** R, int n);
void inicializa(double** mat, double val, int nr, int nc);
void inicializa_vec(double* vec, double val, int n);
double **lee_matriz(char *dir, int *nr, int *nc);
double* crea_vector(int n);
double** crea_matriz(int nr, int nc);
void freeMat(double** mat);

int main(int narg, char **varg)  {
    int nr,nc;

    if(narg != 3) {printf("\nIntroduce la direccion del archivo binario y el numero maximo de operaciones M.\n"); return 0;}
    char *dirmat = varg[1];
    int M = atoi(varg[2]);

    // Leemos los archivos
    double **A = lee_matriz(dirmat,&nr,&nc);
    if(nr != nc) {printf("Error, las matrices en el programa deben ser cuadradas.\n"); freeMat(A); return 0;}
    int n = nr;

    double rem = sqrt(DBL_EPSILON); // tolerancia de raiz de epsilon maquina
    double tolerancia = rem;

    jacobiData data = jacobi(A, M, tolerancia, n);
    printf("El numero de iteraciones hechas fueron: %d\nEl valor de e=|max aij| es: e = %lf\n",data.k,data.e);
    int pvap = 3;
    printf("Los primeros %d eigenvalores de A son:\n",pvap);
    for(int i=0; i<pvap; i++)
        printf("%lf ",data.lm[i]);

    calcula_error_jacobi(A,data.V,data.lm,n);

    free(data.lm);
    freeMat(data.V);
    freeMat(A);
	return 0;
}

void calcula_error_jacobi(double **A, double **V, double *lm, int n)    {
    double **D = crea_matriz(n,n);
    inicializa(D,0,n,n);
    for(int i=0; i<n; i++)
        D[i][i] = lm[i];
    double **X = crea_matriz(n,n); double **Y = crea_matriz(n,n);
    producto_matricial(A,V,X,n);
    producto_matricial(V,D,Y,n);

    double sum=0;
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            sum += pow(X[i][j] - Y[i][j], 2);

    sum = sqrt(sum);
    printf("\nY el error calculado es: ||AV-VD|| = %lf\n", sum);
    freeMat(D); freeMat(X); freeMat(Y);
}

jacobiData jacobi(double **AA, int omax, double tol, int n)  {
	// Hacemos una copia de la matriz porque la vamos a modificar
	double **A = crea_matriz(n,n);
	for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
			A[i][j] = AA[i][j];
    double **V = crea_matriz(n,n);
    int k; double dta;
    indices ind;
    // Iniciamos V = I
    inicializa(V,0,n,n);
    for(int i=0; i<n; i++)
        V[i][i] = 1;

    for(k=1; k<=omax; k++)  {
        // Encontramos indices del mayor valor
        ind = busca_mayor(A,n);
        // Condicion de paro
        if(ind.mx < tol) break;
        // Igualamos el delta a ajj - aii / (2 aij)
        dta = (A[ind.j][ind.j] - A[ind.i][ind.i]) / (2 * A[ind.i][ind.j]);
        // Calculamos t igual al signo(dta) / abs(dta) + raiz de 1 + dta^2
        double signo;
        if(dta < 0) signo = -1;
        else signo = 1;
        //signo = pow(-1,(ind.j - ind.i));
        double t = signo / (abs(dta) + sqrt(1 + dta*dta));
        // Calculamos c igual cosx = 1/ raiz 1+t^2
        double c = 1 / (sqrt(1 + t*t));
        // Calculamos s igual sinx = ct;
        double s = c*t;
        // Creamos la rotacion de givens con los datos c, s, en i y j
        double **G = rotacion_givens(c,s,ind.i,ind.j,n);
        // Calculamos A^k = G(A^k-1)G
        double **X = crea_matriz(n,n);
        producto_matricial(A,G,X,n);
        double **Gt = traspuesta(G,n);
        producto_matricial(Gt,X,A,n);
        // Calculamos V = VG
        producto_matricial(V,G,X,n);
        for(int i=0; i<n; i++)
            for(int j=0; j<n; j++)
                V[i][j] = X[i][j];
        // liberamos memoria auxiliar
        freeMat(X);
        freeMat(Gt);
        freeMat(G);
    }

    // Creamos los datos de jacobi para devolverlos
    jacobiData data;
    data.e = ind.mx;
    data.k = k-1;
    data.V = V;
    data.lm = crea_vector(n);
    for(int i=0; i<n; i++)
        data.lm[i] = A[i][i];
    return data;
}

double** traspuesta(double **A, int n)  {
    double **B = crea_matriz(n,n);
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            B[j][i] = A[i][j];
    return B;
}

double** rotacion_givens(double cose, double sen, int coori, int coorj, int n)    {
    double **G = crea_matriz(n,n);
    inicializa(G,0,n,n);
    for(int i=0; i<n; i++)
        G[i][i] = 1;
    G[coori][coori] = G[coorj][coorj] = cose;
    G[coori][coorj] = sen;
    G[coorj][coori] = -1*sen;
    return G;
}

indices busca_mayor(double **A, int n)  {
    indices ind;
    ind.mx = -1;
    for(int i=0; i<n; i++)
        for(int j=i; j<n; j++)  {
            if(j != i)
                if(ind.mx < abs(A[i][j]))	{
                    ind.mx = abs(A[i][j]);
					ind.i = i;
					ind.j = j;
				}
        }
    return ind;
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
