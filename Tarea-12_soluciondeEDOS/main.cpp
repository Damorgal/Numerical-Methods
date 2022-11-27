#include "EDOS.hpp"

int main(int narg, char **varg)
{
    int n;

    // Leer los datos
    if(narg != 2) {cout << "Introduce el numero n del tamagno de la discretizacion."<< endl; return 0;}
    n = atoi(varg[1]);
		
	// Solucion de PVI descrito abajo
	// Intervalo de solucion
	double a=0, b=5;
	double *y = solve_EDO(a,b,n);
	
	// Hacemos la discretizacion
	double *x = crea_vector(n+1);
	double h = (b-a)/(double)n;
	for(int i=0; i<=n; i++)
		x[i] = a + h*(double)i;
		
	// Error absoluto promedio y error maximo
	double errabs=0, errmax=-1;
	for(int i=0; i<=n; i++)	{
		double aux = abs(y[i]-Y(x[i]));
		errabs += aux;
		if(errmax < aux)
			errmax = aux;
	}
	errabs /= (double)n;
	
	cout<<"El error absoluto promedio es: "<<errabs<<endl;
	cout<<"El error maximo entre las solucion analitica y la numerica es: "<<errmax<<endl;
	cout<<"El valor yn= "<<y[n]<<" y el valor exacto y(xn)= "<< Y(b)<<endl;
	
	free(x);
	free(y);
    return 0;
}

/*** PVI ***/
/*	y'' = -y + xy' - 2xcosx
 * 	y(0) = 0
 * 	y'(0) = 3
La solucion es Y(x) = x + 2sinx

Lo transformamos a:
*  y1' = y2
*  y2' = -y1 + xy2 -2xcosx
*  y1(0) = 0
*  y2(0) = 3

La regla trapezoidal es: yi+1 = yi + h/2*(f(xi,yi)+f(xi+1,yi+1))

El sistema lineal a resolver por metodo del trapecio es:
	| 1     -h/2   |*|y1i+1|   |              y1i + h/2*y2i                 |
	|h/2 1-h/2*xi+1| |y2i+1| = |y2i + h/2*(-y1i+xiy2i-2xicosxi-2xi+1cosxi+1)|
*/
