#include "integracion.hpp"

int main(int narg, char **varg)
{
    int n;
    double a,b,sol_analitica;
    double (*af)(double);
    
    // Leer los datos
    if(narg != 4) {cout << "Introduce los extremos del intervalo (a,b) y el numero n de las divisiones."<< endl; return 0;}
    a = atof(varg[1]);
    b = atof(varg[2]);
    n = atoi(varg[3]);
		
	// Algoritmo
	af = f1;
	double ans = integracion(af,a,b,n);
	cout <<"La integral en ("<<a<<","<<b<<") de f1 es: "<< ans <<endl;
	// La solucion analitica es: sin^2(pi^2)/pi =
	sol_analitica = 0.058937983936327747;
	cout <<"El error relativo de la solucion es: "<<abs(sol_analitica-ans)/sol_analitica <<endl<<endl;
	
	af = f2;
	ans = integracion(af,a,b,n);
	cout <<"La integral en ("<<a<<","<<b<<") de f2 es: "<< ans <<endl;
	// La solucion analitica es: pi-pi^2+pi^4 =
	sol_analitica = 90.68107928650285;
	cout <<"El error relativo de la solucion es: "<<abs(sol_analitica-ans)/sol_analitica <<endl<<endl;

    return 0;
}
