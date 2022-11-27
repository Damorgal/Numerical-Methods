// Funciones hechas para el primer parcial de Metodos numericos
// funciones.cpp hecho para saber lo que hacen las funciones que se utilizan
// Se debe compilar agragando este documento

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

// Funcion que pivotea dos filas de una matriz
void pivoteo_parcial(double** A, int fi, int fj, int nc)	{
	double *aux = crea_vector(nc);
	for(int i=0; i<nc; i++)	{
		aux[i] = A[fi][i];
		A[fi][i] = A[fj][i];
	}
	for(int i=0; i<nc; i++)
		A[fj][i] = aux[i];
	free(aux);
}

// Funcion que pivotea dos filas y columna de una matriz
void pivoteo_total(double** A, int fi, int fj, int ci, int cj, int nr, int nc)	{
	pivoteo_parcial(A,fi,fj,nc);
	double *aux = crea_vector(nc);
	for(int i=0; i<nr; i++)	{
		aux[i] = A[i][ci];
		A[i][ci] = A[i][cj];
	}
	for(int i=0; i<nr; i++)
		A[i][cj] = aux[i];
	free(aux);
}

// Funcion que regresa la matriz P de permutacion para pivoteo parcial o total de matrices cuadradas
double** matriz_pivoteo(int fi, int fj, int n)	{
	double **B = crea_matriz(n,n);
	inicializa_mat(B,0,n,n);
	for(int i=0; i<n; i++)	{
		if(i == fi)
			B[i][fj] = 1;
		else if(i == fj)
			B[i][fi] = 1;
		else
			B[i][i] = 1;
	}
	return B;
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

/**************************** Funciones de errores ****************************/
// Calcula el error de que una matriz por un vector x de la solucion del vector b
void calcula_errorM(double **A, double *b, double *x, int n)  {
    double sum=0;
    for(int i=0; i<n; i++)  {
        double aux=0;
        for(int j=0; j<n; j++)
            aux += A[i][j]*x[j];
        sum += pow(aux - b[i], 2);
    }
    sum = sqrt(sum);
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
            sum += pow(A[i][j] - aux, 2);
        }
    }
    sum = sqrt(sum);
    printf("Y el error calculado de la matriz LU es: %lf\n", sum);
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

void calcula_error_fi(double* fiv, double* x, int n)    {
    double sum=0;
    for(int i=0; i<n; i++)
        sum += (fi(x[i]) - fiv[i]) * (fi(x[i]) - fiv[i]);
    sum = sqrt(sum);
    printf("El error calculado de la solucion analitica con la de ecuaciones es:\n||FI(x) - fi(xi)|| = %lf\n", sum);
}

// Calcula el error de que un a por un vector x de la solucion del vector b
void calcula_error(double *a, double *b, double *x, int n)  {
    double sum=0;
    for(int i=0; i<n; i++)
        sum += pow(a[i]*x[i] - b[i], 2);
    sum = sqrt(sum);
    printf("Y el error calculado de la solucion es: %lf\n", sum);
}

void calcula_error_jacobi(double **A, double **V, double *lm, int n)    {
    double **D = crea_matriz(n,n);
    inicializa_mat(D,0,n,n);
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
    printf("\nY el error calculado es: ||AV-DV|| = %lf\n", sum);
    freeMat(D); freeMat(X); freeMat(Y);
}

/**************************** Funciones de normas ****************************/
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

/****************** Funciones de los metodos: potencia, pot. inversa, jacobi, tridiagonal y valores propios ******************/

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
    inicializa_mat(V,0,n,n);
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
                V[j][i] = X[i][j];
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

double** rotacion_givens(double cos, double sen, int coori, int coorj, int n)    {
    double **G = crea_matriz(n,n);
    inicializa_mat(G,0,n,n);
    for(int i=0; i<n; i++)
        G[i][i] = 1;
    G[coori][coori] = G[coorj][coorj] = cos;
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

double *calcula_d(double *x, double dx, int n)     {
	double k = 0.5;
    double *d = crea_vector(n);
    for(int i=0; i<n; i++)  {
        if(i == 0 || i == n-1) d[i] = fi(x[i]);
        else d[i] = -q(x[i]) * dx * dx / k;
    }
    return d;
}

double* discretizacion(double dx, int n)    {
	double a = 0;
    double *x = crea_vector(n);
    for(int i=0; i<n; i++)
        x[i] = a + i*dx;
    return x;
}

void tiempoTranscurrido(struct timeval t1, struct timeval  t2) {
    double elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;  // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;      // us to ms
    printf("Tiempo transcurrido %6.4f mseg. \n", elapsedTime);
}
