#include "gradiente.hpp"

int main(int narg, char **varg)
{
    int nr,nc,n;
    // Leer los archivos de la matriz A y el vector b
    if(narg != 3) {cout << "Introduce la direccion de la matriz A y el vecctor b."<< endl; return 0;}
    double **A= lee_matriz(varg[1],&nr,&nc);
    if(nr != nc) {cout << "Deben ser matrices cuadradas." << endl; return 0;}
    double *b = lee_arreglo(varg[2],&n);
    if(n != nr) {cout << "La dimension del vector b y la matriz A no coinciden." << endl; return 0;}

    // Metodo del gradiente conjugado
    double tol = sqrt(DBL_EPSILON);
    gradienteData data = gradiente_conjugado(A,b,tol,n);

    // Imprimimos los datos del metodo de gradiente conjugado
    cout << "El tamagno de la matriz es: " << nr <<"x"<< nc <<endl;
    cout << "El numero de iteraciones que se realizaron fue de: " << data.k<< endl;
    cout << "El error que se obtuvo fue de: "<< data.e<< endl;

    double *x = crea_vector(n);
    producto_MatVec(A,data.x,x,n);
    double sum=0;
    for(int i=0; i<n; i++)
        sum += pow(x[i]-b[i],2);
    cout<< "El error ||Ax-b||="<<sqrt(sum)<<endl;

  /*  for(int i=0; i<n; i++)  {
        for(int j=0; j<n; j++)
            cout<<A[i][j]<<" ";
        cout<<endl;
    }
    for(int i=0; i<n; i++)
        cout<<x[i]<<" "<<b[i]<<endl;*/

    free(x);
    // Liberamos memoria
    free(b);
    freeMat(A);
    free(data.x);
    return 0;
}
