#include "funciones.hpp"

int main(int narg, char **varg)
{
	double foco = 0.45;
	place P = inicializa_lugar_vigilancia();
	int ancho = 800, alto = 600; // Pixeles de la camara
	// Debemos cambiar en la funcion la trayectoria deseada
    haz_trayectoria(P,foco,ancho,alto);
	cout<< "Hecho"<<endl;
    return 0;
}
