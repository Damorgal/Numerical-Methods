#include "funciones.hpp"

int main(int narg, char **varg)
{
    int n,M;
    double (*af)(double);
    double x0,tol,h,X;
    
    // numero de divisiones en las integrales
    n = 1000;
    // numero maximo de iteraciones para newton
    M = 1000;
		
	x0 = -2;
	tol = sqrt(DBL_EPSILON);
	h = determina_h(x0,tol);
	cout<<endl<< "Se obtuvo que el h indicado es: h="<<h<<endl;
	cout<< endl<<"Procedamos con newton-rapson para encontrar la raiz"<<endl;
	X = newton_rapson(x0,h,tol,M,n);
	cout<< endl<<"El valor encontrado x* es: "<<X<<" y f(x*) = "<<f(X)<<endl;
	cout<<endl<< "Procedamos a hacer el spline cubico"<<endl;
	haz_spline(X);
	cout<< "Hecho"<<endl;
    return 0;
}
