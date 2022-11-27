// Funciones hechas para la tarea 11 de Integracion Romberg

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

/************************** Funciones para la integgracion de Romberg *****************************/
double f1(double x)	{
	double pi = 3.14159265358979323846264338328;
	return sin(2*pi*x);
}

double f2(double x)	{
	return 4*pow(x,3) - 2*x + 1;
}

double trapecio_recursivo(double (*f)(double), double val, double a, double b, int ind)	{
	double E = pow(2,ind-1), 
	h = (b-a)/pow(2,ind), 
	ans = 0;
	for(int k=1; k<=E; k++)
		ans += f(a + h*(2*k-1));
	ans *= h;
	ans += val/2;
	
	return ans;
}

double integracion(double (*f)(double), double a, double b, int n)	{
	double *R = crea_vector(n+1);
	
	// Calculamos los R(n,0) con regla del trapecio y los guardamos en R
	R[0] = (b-a) * (f(a)+f(b))/2;
	for(int i=1; i<=n; i++)	
		R[i] = trapecio_recursivo(f,R[i-1],a,b,i);
	
	// Extrapolamos con Richardson para obtener R(n,n)
	for(int i=1; i<=n; i++)
		for(int j=i; j<=n; j++)
			R[j-i] = R[j-i+1] + (R[j-i+1]-R[j-i])/(pow(4,j) - 1);
	
	double ans = R[0];
	free(R);
	return ans;
}
