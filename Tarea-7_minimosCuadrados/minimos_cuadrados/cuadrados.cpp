
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

// Inicializa todos las entradas de una matriz a un valor dado
void inicializa_mat(double** mat, double val, int nr, int nc)   {
    for(int i=0; i<nr; i++)
        for(int j=0; j<nc; j++)
            mat[i][j] = val;
}

// Inicializa todos las entradas de una matriz a un valor dado
void inicializa_vec(double* vec, double val, int n)   {
    for(int i=0; i<n; i++)
            vec[i] = val;
}

void inicializa_vector_random(double *vec, int n)  {
    for(int i=0; i<n; i++)
        vec[i] = ((rand()/(double)RAND_MAX)*2)-1;
}

void inicializa_matriz_random(double **A, int n)  {
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            A[i][j] = ((rand()/(double)RAND_MAX)*2)-1;
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

double** traspuesta(double **A, int nr, int nc)  {
    double **B = crea_matriz(nc,nr);
    for(int i=0; i<nr; i++)
        for(int j=0; j<nc; j++)
            B[j][i] = A[i][j];
    return B;
}

// Funcion que hace el producto de una matriz con un vector (P y a) y lo guarda en b (le debes pasar la longitud n de los vectores)
void producto_MatVec(double** P, double* a, double* b, int nr, int nc)   {
    #pragma omp parallel for
	for (int i=0; i<nr; i++) {
        double sum=0;
        //#pragma omp parallel for reduction(+:sum)
        for(int j=0; j<nc; j++)
            sum += P[i][j] * a[j];
        b[i] = sum;
	}
}

// Funcion que hace el producto de dos matrices P y Q y lo guarda en R (le debes pasar la longitud n de los vectores)
void producto_matricial(double** P, double** Q, double** R, int nr, int nc)  {
    #pragma omp parallel for
	for (int i=0; i<nr; i++) {
        //#pragma omp parallel for
        for(int j=0; j<nr; j++)  {
            double sum=0;
            //#pragma omp parallel for reduction(+:sum)
            for(int k=0; k<nc; k++)
                sum += P[i][k] * Q[k][j];
            R[i][j] = sum;
        }
	}
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


/******************************* Funciones para el programa **********************************/
double* crea_valores_xi(int m)  {
    double *x = crea_vector(m);
    double pi = acos(-1);
    for(int i=0; i<m; i++)
        x[i] = i * pi / (m-1);
    return x;
}

double* crea_valores_yi(double* x, int n)  {
    double *y = crea_vector(n);
    for(int i=0; i<n; i++)
        y[i] = sin(x[i]);
    return y;
}

double* minimos_cuadrados(double* x, double* y, int m)   {
    double** A = crea_matriz(m,3);
    // Calculamos FI
    for(int i=0; i<m; i++)  {
        // Funcion f(x) = ax^2 + bx + c, es decir, fi1(x) = x^2, fi2(x) = x, fi3(x) = 1
        A[i][0] = x[i]*x[i];
        A[i][1] = x[i];
        A[i][2] = 1;
    }
    // Calculamos FI^t * FI
    double** At = traspuesta(A,m,3);
    double** B = crea_matriz(3,3);
    producto_matricial(At,A,B,3,m);
    // Calculamos FI^t * y
    double* b = crea_vector(3);
    producto_MatVec(At,y,b,3,m);
    // Resolvemos FI^t * FI * alfa = FI^t * y
    double *alfa = resuelve_ecuacion(B,b,3);
    freeMat(A); freeMat(At);
    freeMat(B); free(b);
    return alfa;
}

double calcula_error_mincuadrados(double* f, double* y, int m)  {
    double sum=0;
    for(int i=0; i<m; i++)
        sum += pow(f[i]-y[i],2);
    sum /= 2;
    return sum;
}
