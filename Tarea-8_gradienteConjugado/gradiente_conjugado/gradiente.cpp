// Funciones hechas para la tarea 8

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

/********* Funciones para metodo del gradiente conjugado *************/
gradienteData gradiente_conjugado(double **A, double *b, double tol, int n) {
    // Iniciamos x=0, r=p=-b y e=sqrt(r*r/n)
    double *x =crea_vector(n), *r = crea_vector(n), *p = crea_vector(n), *w = crea_vector(n);
    inicializa_vec(x,0,n);
    for(int i=0; i<n; i++)	{r[i] = -b[i]; p[i] = -r[i];}
    double e = norma(r,n);
    int k;

    // Algoritmo de gradiente conjugado
    for(k=0; tol<e; k++)  {
        if(k == 2*n ) break;
        producto_MatVec(A,p,w,n);
        double d = producto_punto(r,r,n);
        double alfa = d / producto_punto(p,w,n);
        for(int i=0; i<n; i++)  {
            x[i] += alfa*p[i];
            r[i] += alfa*w[i];
        }
        double beta = producto_punto(r,r,n) / d;
        for(int i=0; i<n; i++)
            p[i] = beta*p[i] - r[i];
        e = norma(r,n);
    }

    // Lenamos datos
    gradienteData data;
    data.k = k;
    data.e = e;
    data.x = x;

    // Liberamos memoria
    free(r); free(p); free(w);
    return data;
}
