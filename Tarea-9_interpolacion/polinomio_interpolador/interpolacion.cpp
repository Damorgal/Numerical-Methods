// Funciones hechas para la tarea 9

// Hecho por Diego Aaron Moreno Galvan

/********************* Funciones generales para la creacion de estructuras ***********************/
// Libera memoria de matriz
void freeMat(double** mat)   {
    free(mat[0]);
    free(mat);
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

    // Lectura del tamagno de la matriz
    fread(nr, sizeof(int), 1, fp);
    fread(nc, sizeof(int), 1, fp);
    // Reservamos memoria
    mat = crea_matriz(*nr, *nc);
    // Lectura de los datos
    fread(mat[0], sizeof(double), (*nr)*(*nc), fp);
    fclose(fp);
    return(mat);
}

void guarda_matriz_txt(double **mat, int nr, int nc, char *cfile) {
    FILE   *fp=fopen(cfile, "wt");
    int     i, j;

    if(!fp) {
        printf("No se puede abrir el archivo\n");
        exit(0);
    }
    printf("Generando el archivo %s ...\n", cfile);

    for(i=0; i<nc; ++i) {
        for(j=0; j<nr; ++j)
            fprintf(fp, "%lf    ", mat[j][i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

/**************************** Funciones de operaciones matriciales y vectoriales ****************************/
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

// Inicializa todos las entradas de una matriz a un valor dado
void inicializa_mat(double** mat, double val, int nr, int nc)   {
    for(int i=0; i<nr; i++)
        for(int j=0; j<nc; j++)
            mat[i][j] = val;
}

// Inicializa todos las entradas de un vector a un valor dado
void inicializa_vec(double* vec, double val, int n)   {
    for(int i=0; i<n; i++)
            vec[i] = val;
}

// Calcula la norma 2 de un vector
double norma(double *vec, int n)   {
    double sum=0;
    for(int i=0; i<n; i++)
        sum += vec[i]*vec[i];
    sum = sqrt(sum);
    return sum;
}

/**************************** Funciones factorizacion LU y solucion de ecuaciones ****************************/
// Resuelve una matriz diagonal (el vector diagonal)
double* resuelve_diagonal(double* a, double* b, int n, int *ind)   {
    double *x = crea_vector(n);
    for(int i=0; i<n; i++)  {
        if(a[i] != 0)
            x[i] = b[i]/a[i];
        else if(b[i] == 0)
            *ind=1;
        else {free(x); return NULL;}
    }
    return x;
}

// Factorizacion por Doolittle, regresa 0 si pudo factorizar, sino 1.
int factoriza_LU(double** A, double** L, double** U, int n)   {
    if(A[0][0] == 0) return 1;
    // Llenamos las matrices L y U de ceros
    inicializa_mat(L,0,n,n);
    inicializa_mat(U,0,n,n);

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

// Lo resolvemos con el caso particular de la gaussiana visto en clase
double* resuelve_tridiagonal(double* d1, double* d2, double* d3, double* d, int n) {
    double *bi = crea_vector(n);
    double *ci = crea_vector(n);
    double *di = crea_vector(n);
    // Procedemos con el algoritmo para crear las entradas modificadas
    for(int i=0; i<n; i++)   {
        if(i == 0)  {
            bi[i] = d2[i];
            ci[i] = d3[i];
            di[i] = d[i];
        }
        else    {
            bi[i] = bi[i-1]*d2[i] - d1[i]*ci[i-1];
            ci[i] = bi[i-1]*d3[i];
            di[i] = bi[i-1]*d[i] - d1[i]*di[i-1];
        }
    }
    // Creamos las n variables fi de la solucion
    double *x = crea_vector(n);
    x[n-1] = d[n-1];
    for(int i=n-2; i>=0; i--)
        x[i] = (di[i] - ci[i]*x[i+1]) / bi[i];

    free(bi);
    free(ci);
    free(di);
    return x;
}

/********* Funciones para la interpolacion de funciones *************/

/* Polinomio interpolador */

double** vandermonde(double *x, int n) {
    double **A = crea_matriz(n+1,n+1);
    for(int i=0; i<=n; i++)
        for(int j=0; j<=n; j++)
            A[i][j] = pow(x[i],j);
    return A;
}

double* coef_pol_interpolador(double *abscisas, double *ordenadas, int n)   {
    double **A = vandermonde(abscisas,n);
    double *c = resuelve_ecuacion(A,ordenadas,n+1);
    freeMat(A);
    return c;
}

double evalua_pol(double *coef, double x, int n)    {
    double ans = 0;
    for(int i=0; i<=n; i++)
        ans += coef[i] * pow(x,i);
    return ans;
}

double** evalua_datos(double *x, double *c, int n)	{
	// Iniciamos constantes y creamos arreglos
	int N = 4*n;
	double dx = (x[n] - x[0])/(double)N;
	double *malla = crea_vector(N+1), *part = crea_vector(N+1);

	// Evaluamos el polinomio
	for(int i=0; i<=N; i++)	{
		part[i] = x[0] + ((double)i * dx);
		malla[i] = evalua_pol(c,part[i],n);
	}

	// Guardamos datos en la matriz
	double **A = crea_matriz(2,N+1);
	for(int j=0; j<N+1; j++)
		A[0][j] = part[j];
	for(int j=0; j<N+1; j++)
		A[1][j] = malla[j];

	// Liberamos memoria
	free(part);
	free(malla);
	return A;
}

/* Polinomio de Lagrange */
double evalua_lagrange(double *x, double val, int ind, int n)    {
    double ans = 1;
    for(int j=0; j<=n; j++)
        if(j != ind)
            ans *= (val-x[j])/(x[ind]-x[j]);
    return ans;
}

double evalua_pol_lagrange(double *x, double *y, double val, int n) {
    double ans = 0;
    for(int j=0; j<=n; j++)
            ans += y[j]*evalua_lagrange(x,val,j,n);
    return ans;
}

double** evalua_datos_lagrage(double *x, double *y, int n)	{
	// Iniciamos constantes y creamos arreglos
	int N = 4*n;
	double dx = (x[n] - x[0])/(double)N;
	double *malla = crea_vector(N+1), *part = crea_vector(N+1);

	// Evaluamos el polinomio
	for(int i=0; i<=N; i++)	{
		part[i] = x[0] + ((double)i * dx);
		malla[i] = evalua_pol_lagrange(x,y,part[i],n);
	}

	// Guardamos datos en la matriz
	double **A = crea_matriz(2,N+1);
	for(int j=0; j<N+1; j++)
		A[0][j] = part[j];
	for(int j=0; j<N+1; j++)
		A[1][j] = malla[j];

	// Liberamos memoria
	free(part);
	free(malla);
	return A;
}

/* Polinomio de diferencias */
double evalua_diferencias(double *x, double val, int ind)   {
    double ans = 1;
    if(ind != 0)
        for(int j=0; j<ind; j++)
            ans *= (val-x[j]);
    return ans;
}

double** matriz_diferencias(double *x, double *y, int n)   {
    double **A = crea_matriz(n+1,n+1);
    for(int i=0; i<n+1; i++)
        A[i][0] = y[i];
    for(int j=1; j<n+1; j++)
        for(int i=0; i<=n-j; i++)
            A[i][j] = (A[i+1][j-1] - A[i][j-1]) / (x[i+j] - x[i]);
    return A;
}

double evalua_pol_diferencias(double *x, double *y, double val, int n)  {
    double ans = 0;
    double **A = matriz_diferencias(x,y,n);
    for(int i=0; i<n+1; i++)
        ans += A[0][i]*evalua_diferencias(x,val,i);
    freeMat(A);
    return ans;
}

double** evalua_datos_diferencias(double *x, double *y, int n)	{
	// Iniciamos constantes y creamos arreglos
	int N = 4*n;
	double dx = (x[n] - x[0])/(double)N;
	double *malla = crea_vector(N+1), *part = crea_vector(N+1);

	// Evaluamos el polinomio
	for(int i=0; i<=N; i++)	{
		part[i] = x[0] + ((double)i * dx);
		malla[i] = evalua_pol_diferencias(x,y,part[i],n);
	}

	// Guardamos datos en la matriz
	double **A = crea_matriz(2,N+1);
	for(int j=0; j<N+1; j++)
		A[0][j] = part[j];
	for(int j=0; j<N+1; j++)
		A[1][j] = malla[j];

	// Liberamos memoria
	free(part);
	free(malla);
	return A;
}
