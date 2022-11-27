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

/************************** Funciones para la integgracion de Romberg *****************************/
double f1(double x)	{
	return exp(-pow(x-5,2)/4);
}

double f(double x)	{
	int n = 8;
    double a = 0, b = x;
    double (*af)(double);
		
	// Algoritmo
	af = f1;
	double ans = integracion(af,a,b,n);
	//cout <<"La integral en ("<<a<<","<<b<<") de f1 es: "<< ans <<endl;
	
	return pow(x-3,3)/20 + 10*ans - 3;
}

double derf(double x, double h)	{
	return (f(x+h) - f(x-h))/(2*h);
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

double newton_rapson(double x, double h, double tol, int M, int n)	{
	double ans;
	for(int i=0; i<M; i++)	{
		ans = f(x);
		if(abs(ans) < tol)
			break;
		cout<< "Indice k = "<< i<<endl;
		cout<< "Valor de la funcion f(xk) = "<< ans<< endl;
		cout<< "El punto xk = "<< x<< endl;
		x = x - (f(x)/derf(x,h));
	}
	return x;
}

double determina_h(double x, double tol)	{
	double ans,h1,h0;
	for(int i=0; i<100; i++)	{
		h0 = 0.01-(double)i/10000;
		h1 = 0.01-((double)i+1)/10000;
		ans = derf(x,h1) - derf(x,h0);
		if(abs(ans) < tol)
			break;
	}
	return h0;
}

void haz_spline(double x)	{
	int m = 60;
	// Hacemos la particion
	double h = 8.0/10.0;
	double *px = crea_vector(11);
	double *py = crea_vector(11);
	double **A = crea_matriz(2,11);
	for(int i=0; i<=10; i++)	{
		px[i] = x-4 + i*h;
		A[0][i] = px[i];
	}
	for(int i=0; i<=10; i++) {
		py[i] = f(px[i]);
		A[1][i] = py[i];
	}
	
	// Algoritmo
    double *M = coef_spline_cubico(px,py,10);
	double **datos = evalua_datos_spline(px,py,M,10,m);
	
	// Imprimimos la matriz de datos en txt para graficar
	char s1[] = "puntos.txt";
	guarda_matriz_txt(A,2,11,s1);
	char s[] = "puntos_splines.txt";
	guarda_matriz_txt(datos,2,m+1,s);

    // Liberamos memoria
    freeMat(datos); freeMat(A);
    free(M);
    free(px); free(py);
    return;
}

/********************************** Splines cubicos **********************************/
double* coef_spline_cubico(double *px, double *py, int n)	{
	// Creamos las diagonales para la matriz tridiagonal
	double *d1 =crea_vector(n+1), *d2 = crea_vector(n+1), *d3 = crea_vector(n+1);
	inicializa_vec(d1,0,n+1);
	inicializa_vec(d2,2,n+1);
	inicializa_vec(d3,0,n+1);
	double *ordenadas = crea_vector(n+1);
	inicializa_vec(ordenadas,0,n+1);
	
	// Creamos sistema de ecuaciones y ordenadas
	double h1,h2,mu,lm,d;
	for(int i=1; i<n; i++)	{
		h1 = px[i] - px[i-1];
		h2 = px[i+1] - px[i];
		mu = h1/(h1+h2);
		lm = h2/(h1+h2);
		d = 6/(h1+h2) * ((py[i+1]-py[i])/h2 - (py[i]-py[i-1])/h1);
		d1[i] = mu;
		d3[i] = lm;
		ordenadas[i] = d;
	}
	
	// Resolvemos el sistema tridiagonal
	double *x = resuelve_tridiagonal(d1,d2,d3,ordenadas,n+1);
	
	// Liberamos memoria
	free(ordenadas);
	free(d1); free(d2); free(d3);
	
	return x;
}

double evalua_pol_cubico(double *px, double *py, double *coef, double val, int n)	{
	// Encontramos intervalo
	int ind;
	for(int i=0; i<n; i++)
		if(val >= px[i] && val <= px[i+1])	{
			ind = i;
			break;
		}
	
	// Calculamos preliminares
	double h,CC,C;
	h = px[ind+1] - px[ind];
	CC = py[ind] - coef[ind]*h*h/6;
	C = (py[ind+1]-py[ind])/h - h/6*(coef[ind+1]-coef[ind]);
	
	// Evaluamos polinomio
	double ans = coef[ind]*(pow(px[ind+1]-val,3))/(6*h) + 
	coef[ind+1]*pow(val-px[ind],3)/(6*h) + 
	C*(val-px[ind]) + CC;
	
	return ans;
}

double** evalua_datos_spline(double *px, double *py, double *coef, int n, int m)	{
	double **A =crea_matriz(2,m+1);
	double h = (px[n]-px[0])/(double)m;
	for(int i=0; i<m+1; i++) {
		A[0][i] = px[0]+(double)i*h;
		A[1][i] = evalua_pol_cubico(px,py,coef,A[0][i],n);
	}	
	return A;
}
