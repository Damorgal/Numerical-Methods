/*  Metodos numericos
    Tarea 1
    Diego Aaron Moreno Galvan
*/

/*  Para compilar:
    gcc -o Tarea1 main.cpp  -lgsl -lgslcblas -lm

    Para ejecutar:
    ./Tarea1
*/

#include <stdio.h>
#include <float.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_ieee_utils.h>
#include <fenv.h>

void ValoresEpsilon();
void FuncionRaiz(double x);

int main(void)  {

    ValoresEpsilon();
    FuncionRaiz(2);
    FuncionRaiz(100000000000);

    return 0;
}

void ValoresEpsilon()   {
    double epsilon, aux;

    // Calculo dado por DBL_EPSILON
    epsilon = DBL_EPSILON;
    printf("\nEpsilon dado por DBL_EPSILON: \nDecimal = ");
    printf("%28.20e\nBinario = ",epsilon);
    gsl_ieee_printf_double(& epsilon);

    // Resultado de suma de 1.0 + epsilon
    printf("\n\nSuma de 1.0 + DBL_EPSILON: \nDecimal = ");
    aux = 1.0 + epsilon;
    printf("%28.20e\nBinario = ",aux);
    gsl_ieee_printf_double(& aux);

    // Resultado de 1.0 - epsilon/2
    printf("\n\nResta de 1.0 - DBL_EPSILON / 2: \nDecimal = ");
    aux = 1.0 - (epsilon/2.0);
    printf("%28.20e\nBinario = ",aux);
    gsl_ieee_printf_double(& aux);
}

void FuncionRaiz(double x)  {
    double x0 = x;

    // Imprimir valores en Decimal y Binario
    printf("\n\nLa tabla de valores para x = %g:\n  Valor             Decimal                                      Binario\n",x);
    for(int i=1; i<=60; ++i)
    {
        // Imprimir Valor
        if(i<10)
            printf("    %d    ",i);
        else
            printf("   %d    ",i);

        // Imprimir Decimal y binario
        printf("%28.20e   ",x);
        gsl_ieee_printf_double(& x); printf("\n");

        x = sqrt(x);
        if(x==x0) break;
        else x0 = x;
    }
}
