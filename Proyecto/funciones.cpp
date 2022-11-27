// Funciones hechas para el proyecto final

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

double evalua_pol_cubico(double *px, double *py, double *coef, double val, int n, int ind)	{
	/*// Encontramos intervalo
	int ind;
	for(int i=0; i<n; i++)
		if(val >= px[i] && val <= px[i+1])	{
			ind = i;
			break;
		}
	*/
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
/*
double** evalua_datos_spline(double *px, double *py, double *coef, int n, int m)	{
	double **A =crea_matriz(2,m+1);
	double h = (px[n]-px[0])/(double)m;
	for(int i=0; i<m+1; i++) {
		A[0][i] = px[0]+(double)i*h;
		A[1][i] = evalua_pol_cubico(px,py,coef,A[0][i],n);
	}	
	return A;
}*/

/* Polinomio de Lagrange */
double evalua_lagrange(double *x, double val, int ind, int n)    {
    double ans = 1.0;
    for(int j=0; j<=n; j++)
        if(j != ind)
            ans *= (val-x[j])/(x[ind]-x[j]);
    return ans;
}

double evalua_pol_lagrange(double *x, double *y, double val, int n) {
    double ans = 0.0;
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

/****************************************************** Funciones para proyecto ***********************************************************/
place inicializa_lugar_vigilancia()	{
	place P;
	P.largo = 10.0;
	P.ancho = 10.0;
	P.poscam1.x = 0.0; P.poscam2.x = P.ancho;
	P.poscam1.y = P.poscam2.y = 0.0;
	P.poscam1.z = P.poscam2.z = 3.0;
	P.apcam1.x = P.apcam2.x = P.ancho/2.0;
	P.apcam1.y = P.apcam2.y = P.largo/2.0;
	P.apcam1.z = P.apcam2.z = 0.0;
	return P;
}

void normaliza(punto* vec)	{
	double norma = vec->x*vec->x + vec->y*vec->y + vec->z*vec->z;
	norma = sqrt(norma);
	vec->x /= norma; vec->y /= norma; vec->z /= norma;
}

double norma(punto *vec)	{
	double ans = vec->x*vec->x + vec->y*vec->y + vec->z*vec->z;
	ans = sqrt(ans);
	return ans;
}
/*
punto trayectoria(double t)	{
	punto ans;
	ans.x = sin(t)+4;
	ans.y = 2*t + 2;
	ans.z = sin(t) + 1;
	return ans;
}

punto trayectoria(double t)	{
	punto ans;
	ans.x = t;
	ans.y = 2*t + 2;
	ans.z = sin(t) + 1;
	return ans;
}
*/
punto trayectoria(double t)	{
	punto ans;
	ans.x = t;
	ans.y = cos(t)+2;
	ans.z = sin(2*t) + 1;
	return ans;
}

double angulo_vectorial(punto* u, punto* v)	{
	double nu = norma(u), nv = norma(v);
	double pp;
	pp = u->x*v->x + u->y*v->y + u->z*v->z;
	return acos(pp/(nu*nv));
}

double* proyeccion_cam1(punto pos, place P, double foco, int ancho, int alto, double h)	{
	double factor;
	punto vecap,vec,pto;
	vecap.x = P.poscam1.x - P.apcam1.x;
	vecap.y = P.poscam1.y - P.apcam1.y;
	vecap.z = P.poscam1.z - P.apcam1.z;
	factor = foco / norma(&vecap);
	vec.x = P.poscam1.x - pos.x;
	vec.y = P.poscam1.y - pos.y;
	vec.z = P.poscam1.z - pos.z;
	//normaliza(&vecap);
	//normaliza(&vec);
	pto.x = -P.apcam1.x*factor; pto.x = P.poscam1.x - pto.x;
	pto.y = -P.apcam1.y*factor; pto.y = P.poscam1.y - pto.y;
	pto.z = P.poscam1.z*factor; pto.z = P.poscam1.z - pto.z;
	double **A = crea_matriz(3,3); inicializa_mat(A,0,3,3);
	double *b = crea_vector(3); inicializa_vec(b,0,3);
	A[0][0] = vecap.x; A[0][1] = vecap.y; A[0][2] = vecap.z;
	A[1][0] = 1/vec.x; A[1][1] = -1/vec.y;
	A[2][1] = 1/vec.y; A[2][2] = -1/vec.z;
	b[0] = vecap.x*pto.x + vecap.y*pto.y + vecap.z*pto.z;
	b[1] = pos.x/vec.x - pos.y/vec.y;
	b[2] = pos.y/vec.y - pos.z/vec.z;
	double *proyeccion = resuelve_ecuacion(A,b,3);
	punto ans;
	ans.x = proyeccion[0]; ans.y = proyeccion[1]; ans.z = proyeccion[2];
	
	A[0][0] = 1.0; A[0][1] = 0.0; A[0][2] = 0.0;
	A[1][0] = 0.0; A[1][1] = 1.0; A[1][2] = 0.0; 
	A[2][0] = vecap.x; A[2][1] = vecap.y; A[2][2] = vecap.z;
	b[0] = 0.0;
	b[1] = 0.0;
	b[2] = vecap.x*pto.x + vecap.y*pto.y + vecap.z*pto.z;
	double *aux = resuelve_ecuacion(A,b,3);
	punto eje;
	eje.x = aux[0]; eje.y = aux[1]; eje.z = aux[2];
	
	eje.x -= pto.x; eje.y -= pto.y; eje.z -= pto.z;
	ans.x -= pto.x; ans.y -= pto.y; ans.z -= pto.z;
	
	double mag1 = norma(&eje), mag2 = norma(&ans);
	double angle = angulo_vectorial(&eje,&ans);
	punto t,q; t.x = eje.x; t.y = eje.y; t.z = 0.0;
	q.x = 0.0; q.y = -1.0; q.z = 0.0;
	double an = angulo_vectorial(&t,&q);
	t.x = ans.x; t.y = ans.y;
	double q1 = t.x*cos(an) - t.y*sin(an);
	if(q1 < 0)
		angle *= -1;
	eje.x = 0.0; eje.y = -mag1; eje.z = 0.0;
	ans.x = 0.0; ans.y = -mag2; ans.z = 0.0;
	double p1 = ans.x*cos(angle) - ans.y*sin(angle);
	double p2 = ans.x*sin(angle) + ans.y*cos(angle);
	ans.x = p1; ans.y = p2;
	
	// Intervalos de los pixeles de camara con el h dado
	double *x = crea_vector(ancho+1);
	double *y = crea_vector(alto+1);
	for(int i=0; i<=ancho; i++)
		x[i] = i*h;
	for(int i=0; i<=alto; i++)
		y[i] = i*h;
		
	t.x = x[ancho/2];
	t.y = y[alto/2];
	ans.x += t.x; ans.y += t.y;
	
	// Checamos en que pixel cayo
	double *answ = crea_vector(2);
	for(int i=0; i<ancho; i++)
		if(ans.x >= x[i] && ans.x <= x[i+1])	{
			answ[0] = (double)i-h/2.0;
			break;
		}
	for(int i=0; i<alto; i++)
		if(ans.y >= y[i] && ans.y <= y[i+1])	{
			answ[1] = (double)i-h/2.0;
			break;
		}
	
	// Borrar memoria
	freeMat(A);
	free(aux); free(b); free(x);
	free(y); free(proyeccion);
	return answ;
}

double* proyeccion_cam2(punto pos, place P, double foco, int ancho, int alto, double h)	{
	double factor;

	punto vecap,vec,pto;
	vecap.x = P.poscam2.x - P.apcam2.x;
	vecap.y = P.poscam2.y - P.apcam2.y;
	vecap.z = P.poscam2.z - P.apcam2.z;
	factor = foco / norma(&vecap);
	vec.x = P.poscam2.x - pos.x;
	vec.y = P.poscam2.y - pos.y;
	vec.z = P.poscam2.z - pos.z;
	//normaliza(&vecap);
	//normaliza(&vec);
	pto.x = P.apcam2.x*factor; pto.x = P.poscam2.x - pto.x;
	pto.y = -P.apcam2.y*factor; pto.y = P.poscam2.y - pto.y;
	pto.z = P.poscam2.z*factor; pto.z = P.poscam2.z - pto.z;
	double **A = crea_matriz(3,3); inicializa_mat(A,0,3,3);
	double *b = crea_vector(3); inicializa_vec(b,0,3);
	A[0][0] = vecap.x; A[0][1] = vecap.y; A[0][2] = vecap.z;
	A[1][0] = 1/vec.x; A[1][1] = -1/vec.y;
	A[2][1] = 1/vec.y; A[2][2] = -1/vec.z;
	b[0] = vecap.x*pto.x + vecap.y*pto.y + vecap.z*pto.z;
	b[1] = pos.x/vec.x - pos.y/vec.y;
	b[2] = pos.y/vec.y - pos.z/vec.z;
	double *proyeccion = resuelve_ecuacion(A,b,3);
	punto ans;
	ans.x = proyeccion[0]; ans.y = proyeccion[1]; ans.z = proyeccion[2];
	
	A[0][0] = 1.0; A[0][1] = 0.0; A[0][2] = 0.0;
	A[1][0] = 0.0; A[1][1] = 1.0; A[1][2] = 0.0; 
	A[2][0] = vecap.x; A[2][1] = vecap.y; A[2][2] = vecap.z;
	b[0] = P.poscam2.x;
	b[1] = 0.0;
	b[2] = vecap.x*pto.x + vecap.y*pto.y + vecap.z*pto.z;
	double *aux = resuelve_ecuacion(A,b,3);
	punto eje;
	eje.x = aux[0]; eje.y = aux[1]; eje.z = aux[2];
	
	eje.x -= pto.x; eje.y -= pto.y; eje.z -= pto.z;
	ans.x -= pto.x; ans.y -= pto.y; ans.z -= pto.z;
	
	double mag1 = norma(&eje), mag2 = norma(&ans);
	double angle = angulo_vectorial(&eje,&ans);
	punto t,q; t.x = eje.x; t.y = eje.y; t.z = 0.0;
	q.x = 0.0; q.y = -1.0; q.z = 0.0;
	double an = -angulo_vectorial(&t,&q);
	t.x = ans.x; t.y = ans.y;
	double q1 = t.x*cos(an) - t.y*sin(an);
	if(q1 < 0)
		angle *= -1;
	eje.x = 0.0; eje.y = -mag1; eje.z = 0.0;
	ans.x = 0.0; ans.y = -mag2; ans.z = 0.0;
	double p1 = ans.x*cos(angle) - ans.y*sin(angle);
	double p2 = ans.x*sin(angle) + ans.y*cos(angle);
	ans.x = p1; ans.y = p2;
	
	// Intervalos de los pixeles de camara con el h dado
	double *x = crea_vector(ancho+1);
	double *y = crea_vector(alto+1);
	for(int i=0; i<=ancho; i++)
		x[i] = i*h;
	for(int i=0; i<=alto; i++)
		y[i] = i*h;
		
	t.x = x[ancho/2];
	t.y = y[alto/2] - mag1 - eje.y;
	ans.x += t.x; ans.y += t.y;
	
	// Checamos en que pixel cayo
	double *answ = crea_vector(2);
	for(int i=0; i<ancho; i++)
		if(ans.x >= x[i] && ans.x <= x[i+1])	{
			answ[0] = (double)i-h/2.0;
			break;
		}
	for(int i=0; i<alto; i++)
		if(ans.y >= y[i] && ans.y <= y[i+1])	{
			answ[1] = (double)i-h/2.0;
			break;
		}
	
	// Borrar memoria
	freeMat(A);
	free(aux); free(b); free(x);
	free(y); free(proyeccion);
	return answ;
}

double error_extrapolacion(double x, double y, double z)	{
	punto p = trayectoria(x);
	double ans = pow(p.y-y,2) + pow(p.z-z,2);
	return sqrt(ans);
}

void extrapola_splines(double **T, int m)	{
	double delta = .2,x,y,z;
	x = T[0][m-1] + delta;
	double *u = crea_vector(m), *v = crea_vector(m);
	
	for(int i=0; i<m; i++)	{
		u[i] = T[0][i];
		v[i] = T[1][i];
	}
	double *M = coef_spline_cubico(u,v,m-1);
	y = evalua_pol_cubico(u,v,M,x,m-1,m-2); free(M);
	
	for(int i=0; i<m; i++)	{
		u[i] = T[0][i];
		v[i] = T[2][i];
	}
	double *N = coef_spline_cubico(u,v,m-1);
	z = evalua_pol_cubico(u,v,N,x,m-1,m-2); free(N);
	double **TE = crea_matriz(3,m+1);
	TE[0][m] = x; TE[1][m] = y; TE[2][m] = z;
	for(int i=0; i<m; i++)	{
		TE[0][i] = T[0][i];
		TE[1][i] = T[1][i];
		TE[2][i] = T[2][i];
	} 
	char s[] = "extrapolacion_splines.txt";
	guarda_matriz_txt(TE,3,m+1,s);
	freeMat(TE);
	
	cout<<"El error de la extrapolacion con splines es: "<<error_extrapolacion(x,y,z)<<endl;
}

void extrapola(double **T, int m)	{
	double delta = .5,x,y,z; int n;
	n = 0;
	x = T[0][m-1] + delta;
	
	// Extrapolacion lineal
	y = T[1][m-2] + (x-T[0][m-2])*(T[1][m-1]-T[1][m-2])/(T[0][m-1]-T[0][m-2]);
	z = T[2][m-2] + (x-T[0][m-2])*(T[2][m-1]-T[2][m-2])/(T[0][m-1]-T[0][m-2]);
	
	// Extrapolacion con polinomios de Lagrange
	/**double *u = crea_vector(m-n), *v = crea_vector(m-n);
	for(int i=n; i<m; i++)	{
		u[i-n] = T[0][i];
		v[i-n] = T[1][i];
	}
	y = evalua_pol_lagrange(u,v,x,m-n-1);
	for(int i=n; i<m; i++)	{
		u[i-n] = T[0][i];
		v[i-n] = T[2][i];
	}
	z = evalua_pol_lagrange(u,v,x,m-n-1);**/
	
	// Extrapolacion de richardson solo para los ultimos 2 puntos
	/**double *R = crea_vector(2);
	for(int i=0; i<2; i++)
		R[i] = T[1][m-2+i];
	
	// Extrapolamos con Richardson 
	R[0] = R[1] + (R[1]-R[0])/(3);
	y = R[0];
	for(int i=0; i<2; i++)
		R[i] = T[2][m-2+i];
	
	// Extrapolamos con Richardson 
	R[0] = R[1] + (R[1]-R[0])/(3);
	z = R[0];
	free(R);**/
	
	// Extrapolacion con polinomios de Newton
	/**double *u = crea_vector(m-n), *v = crea_vector(m-n);
	for(int i=n; i<m; i++)	{
		u[i-n] = T[0][i];
		v[i-n] = T[1][i];
	}
	y = evalua_pol_diferencias(u,v,x,m-n-1);
	for(int i=n; i<m; i++)	{
		u[i-n] = T[0][i];
		v[i-n] = T[2][i];
	}
	z = evalua_pol_diferencias(u,v,x,m-n-1);**/
	
	double **TE = crea_matriz(3,m+1);
	TE[0][m] = x; TE[1][m] = y; TE[2][m] = z;
	for(int i=0; i<m; i++)	{
		TE[0][i] = T[0][i];
		TE[1][i] = T[1][i];
		TE[2][i] = T[2][i];
	} 
	char s[] = "extrapolacion_lineal.txt";
	guarda_matriz_txt(TE,3,m+1,s);
	freeMat(TE);
	
	cout<<"El error de la extrapolacion lineal es: "<<error_extrapolacion(x,y,z)<<endl;
}

double error_trayectoria(double **T, int m)	{
	double ans=0.0;
	double **t = crea_matriz(3,m);
	FILE *f = fopen("trayectoria.txt", "rt");
	for(int i=0; i<m; i++)	{
		fscanf(f,"%lf ",&t[0][i]);
		fscanf(f,"%lf ",&t[1][i]);
		fscanf(f,"%lf ",&t[2][i]);
	} fclose(f);
	for(int i=0; i<m; i++)	{
		ans += pow(t[0][i]-T[0][i],2);
		ans += pow(t[1][i]-T[1][i],2);
		ans += pow(t[2][i]-T[2][i],2);
	} ans = sqrt(ans);
	return ans;
}

void reconstruye_trayectoria(double** C1, double** C2, place P, double m)	{
	punto v,u,t;
	double **T = crea_matriz(3,m),error;
	// Por la izquierda
	for(int i=0; i<m; i++)	{
		u.x = C1[0][i] - P.poscam1.x; v.x = C2[0][i] - P.poscam2.x;
		u.y = C1[1][i] - P.poscam1.y; v.y = C2[1][i] - P.poscam2.y;
		u.z = C1[2][i] - P.poscam1.z; v.z = C2[2][i] - P.poscam2.z;
		t.x = -u.y; t.y = u.x; t.z = 0.0;
		double **A = crea_matriz(3,3); inicializa_mat(A,0,3,3);
		double *b = crea_vector(3); inicializa_vec(b,0,3);
		
		A[0][0] = t.x; A[0][1] = t.y;
		b[0] = t.x*u.x + t.y*u.y;
		A[1][0] = 1/v.x; A[1][1] = -1/v.y;
		b[1] = C2[0][i]/v.x - C2[1][i]/v.y;
		A[2][0] = 1/v.x; A[2][2] = -1/v.z;
		b[2] = C2[0][i]/v.x - C2[2][i]/v.z;
		double *x = resuelve_ecuacion(A,b,3);
		T[0][i] = x[0];
		T[1][i] = x[1];
		T[2][i] = x[2];
		free(x); free(b); 
		freeMat(A); 
	}
	char s1[] = "trayectoria_reconstruida1.txt";
	guarda_matriz_txt(T,3,m,s1);
	cout<<"El error de reconstruccion por la izquierda es: "<<error_trayectoria(T,m)<<endl;

	double **T2 = crea_matriz(3,m);
	// Por la derecha
	for(int i=0; i<m; i++)	{
		u.x = C1[0][i] - P.poscam1.x; v.x = C2[0][i] - P.poscam2.x;
		u.y = C1[1][i] - P.poscam1.y; v.y = C2[1][i] - P.poscam2.y;
		u.z = C1[2][i] - P.poscam1.z; v.z = C2[2][i] - P.poscam2.z;
		t.x = -v.y; t.y = v.x; t.z = 0.0;
		double **A = crea_matriz(3,3); inicializa_mat(A,0,3,3);
		double *b = crea_vector(3); inicializa_vec(b,0,3);
		
		A[0][0] = t.x; A[0][1] = t.y;
		b[0] = t.x*C2[0][i] + t.y*C2[1][i];
		A[1][0] = 1/u.x; A[1][1] = -1/u.y;
		b[1] = C1[0][i]/u.x - C1[1][i]/u.y;
		A[2][0] = 1/u.x; A[2][2] = -1/u.z;
		b[2] = C1[0][i]/u.x - C1[2][i]/u.z;
		double *x = resuelve_ecuacion(A,b,3);
		T2[0][i] = x[0];
		T2[1][i] = x[1];
		T2[2][i] = x[2];
		free(x); free(b); 
		freeMat(A); 
	}
	char s2[] = "trayectoria_reconstruida2.txt";
	guarda_matriz_txt(T2,3,m,s2);
	cout<<"El error de reconstruccion por la derecha es: "<<error_trayectoria(T2,m)<<endl;
	
	// Hacemos el promedio de las dos trayectoria
	for(int i=0; i<m; i++)	{
		T[0][i] = (T[0][i] + T2[0][i])/2.0;
		T[1][i] = (T[1][i] + T2[1][i])/2.0;
		T[2][i] = (T[2][i] + T2[2][i])/2.0;
	}
	char s[] = "trayectoria_reconstruida.txt";
	guarda_matriz_txt(T,3,m,s);
	cout<<"El error de reconstruccion promedio: "<<error_trayectoria(T,m)<<endl;
	
	extrapola(T,m);
	extrapola_splines(T,m);
	freeMat(T); freeMat(T2);
}

void proyecta_imagenes(place P, double foco, int ancho, int alto, double h, double** C1, double** C2, int m)	{
	punto u,v,t,s;
	double pi = 3.14159265,factor,angle,p1,p2;
	// Intervalos de los pixeles de camara con el h dado
	double *x = crea_vector(ancho+1);
	double *y = crea_vector(alto+1);
	for(int i=0; i<=ancho; i++)
		x[i] = i*h;
	for(int i=0; i<=alto; i++)
		y[i] = i*h;
	
	// Checamos el punto correspondiente del pixel
	for(int i=0; i<m; i++)	{
		int f = floor(C1[0][i]);
		C1[0][i] = (x[f] + x[f+1])/2.0;
		f = floor(C2[0][i]);
		C2[0][i] = (x[f] + x[f+1])/2.0;
		f = floor(C1[1][i]);
		C1[1][i] = (y[f+1] + y[f])/2.0;
		f = floor(C2[1][i]);
		C2[1][i] = (y[f+1] + y[f])/2.0;
	}
	char s2[] = "proyeccion_cam1.txt";
	guarda_matriz_txt(C1,2,m,s2);
	
	char s3[] = "proyeccion_cam2.txt";
	guarda_matriz_txt(C2,2,m,s3);
	
	// Trasladamos puntos a camara
	double **C13d = crea_matriz(3,m);
	double **C23d = crea_matriz(3,m);
	for(int i=0; i<m; i++)	{
		// Para la camara 1
		u.x = C1[0][i] - x[ancho/2];
		u.y = C1[1][i] - y[alto/2];   u.z=0.0;
		t.x = P.apcam1.x - P.poscam1.x;
		t.y = P.apcam1.y - P.poscam1.y;
		t.z = P.apcam1.z - P.poscam1.z;
		angle = pi/2 - angulo_vectorial(&t,&P.apcam1);
		p1 = u.y*cos(angle) - u.z*sin(angle);
		p2 = u.y*sin(angle) + u.z*cos(angle);
		u.y = p1; u.z = p2;
		t.x = 0.0; t.y = 1.0; t.z = 0.0;
		angle = -angulo_vectorial(&P.apcam1,&t);
		p1 = u.x*cos(angle) - u.y*sin(angle);
		p2 = u.x*sin(angle) + u.y*cos(angle);
		u.x = p1; u.y = p2;
		
		t.x = P.poscam1.x - P.apcam1.x;
		t.y = P.poscam1.y - P.apcam1.y;
		t.z = P.poscam1.z - P.apcam1.z;
		factor = foco / norma(&t);
		s.x = -P.apcam1.x*factor; s.x = P.poscam1.x - s.x;
		s.y = -P.apcam1.y*factor; s.y = P.poscam1.y - s.y;
		s.z = P.poscam1.z*factor; s.z = P.poscam1.z - s.z;
		v.x = s.x + u.x; v.y = s.y + u.y; v.z = s.z + u.z;  // v ya es el punto de la imagen en la camara
		C13d[0][i] = v.x; C13d[1][i] = v.y; C13d[2][i] = v.z; // Los guardamos en la matriz por mientras
		
		// Para la camara 2
		u.x = C2[0][i] - x[ancho/2];
		u.y = C2[1][i] - y[alto/2];   u.z=0.0;
		t.x = P.apcam2.x - P.poscam2.x;
		t.y = P.apcam2.y - P.poscam2.y;
		t.z = P.apcam2.z - P.poscam2.z;
		s.x = t.x; s.y = t.y; s.z=0.0;
		angle = pi/2 - angulo_vectorial(&t,&s);
		p1 = u.y*cos(angle) - u.z*sin(angle);
		p2 = u.y*sin(angle) + u.z*cos(angle);
		u.y = p1; u.z = p2;
		t.x = 0.0; t.y = 1.0; t.z = 0.0;
		angle = angulo_vectorial(&P.apcam2,&t);
		p1 = u.x*cos(angle) - u.y*sin(angle);
		p2 = u.x*sin(angle) + u.y*cos(angle);
		u.x = p1; u.y = p2;
		
		t.x = P.poscam2.x - P.apcam2.x;
		t.y = P.poscam2.y - P.apcam2.y;
		t.z = P.poscam2.z - P.apcam2.z;
		factor = foco / norma(&t);
		s.x = P.apcam2.x*factor; s.x = P.poscam2.x - s.x;
		s.y = -P.apcam2.y*factor; s.y = P.poscam2.y - s.y;
		s.z = P.poscam2.z*factor; s.z = P.poscam2.z - s.z;
		v.x = s.x + u.x; v.y = s.y + u.y; v.z = s.z + u.z;  // v ya es el punto de la imagen en la camara
		C23d[0][i] = v.x; C23d[1][i] = v.y; C23d[2][i] = v.z; // Los guardamos en la matriz por mientras
	}
	reconstruye_trayectoria(C13d,C23d,P,m);
	char s4[] = "proyeccion_real_cam1.txt";
	guarda_matriz_txt(C13d,3,m,s4);
	char s5[] = "proyeccion_real_cam2.txt";
	guarda_matriz_txt(C23d,3,m,s5);
	
	freeMat(C1); freeMat(C13d);
	freeMat(C2); freeMat(C23d);
}

void haz_trayectoria(place P, double foco, int ancho, int alto)	{
	double h1 = 0.4, h2 = 0.001;
	int m = 20;
	punto pto; 
	double **A = crea_matriz(3,m);
	double **B = crea_matriz(2,m);
	double **C = crea_matriz(2,m);
	for(int i=0; i<m; i++)	{
		pto = trayectoria(((double)i+0.001)*h1);
		double *ans1 = proyeccion_cam1(pto,P,foco,ancho,alto,h2);
		double *ans2 = proyeccion_cam2(pto,P,foco,ancho,alto,h2);
		A[0][i] = pto.x;
		A[1][i] = pto.y;
		A[2][i] = pto.z;
		B[0][i] = ans1[0];
		B[1][i] = ans1[1];
		C[0][i] = ans2[0];
		C[1][i] = ans2[1];
		free(ans1); free(ans2);
	}
	char s[] = "trayectoria.txt";
	guarda_matriz_txt(A,3,m,s);
	freeMat(A);
	
	char s2[] = "imagen_cam1.txt";
	guarda_matriz_txt(B,2,m,s2);
	
	char s3[] = "imagen_cam2.txt";
	guarda_matriz_txt(C,2,m,s3);
	
	proyecta_imagenes(P,foco,ancho,alto,h2,B,C,m);
}
