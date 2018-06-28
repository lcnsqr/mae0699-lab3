#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define RAND ((double)(rand() >> 1)/((RAND_MAX >> 1) + 1))

// Função objetivo para encontrar a raiz
double f(double a, double b, double x){
	return (1-a*x)/(1-b*x) - exp((b-a)*x);
}

// Derivada da função objetivo
double df(double a, double b, double x){
	return (b-a)*(pow(1-b*x,-2)-exp((b-a)*x));
}

// Determinar o valor de k 
double k(double a, double b, double x){
	return x/(exp(b*x)-exp(a*x));
}

int main(int argc, char** argv){
	srand((unsigned int)time(NULL));

	// Constantes da função
	double a, b;
	a = atof(argv[1]);
	b = atof(argv[2]);
	// Tolerância
	double T;
	T = atof(argv[3]);
	// Número máximo de iterações
	int N;
	N = atoi(argv[4]);
	// Solução (lambda)
	double p;
	// Contagem de iterações
	int i;

	// Aproximação inicial
	double p0;
	// Iterar até encontrar um lambda não muito próximo de zero
	while (1){
		// Variar a aproximação inicial em [-10,10]
		p0 = -50.0 + 100.0 * RAND;
		printf("Começando em %.6f\n", p0);
		i = 1;
		while ( i < N ){
			p = p0 - f(a, b, p0)/df(a, b, p0);
			//printf("x = %.36f\n", p);
			if ( fabs(p-p0) < T && p > 0.00001 ){
				printf("l = %.36f\n", p);
				//printf("f(x) = %.36f\n", f(a, b, p));
				//printf("Iterações: %d\n", i);
				printf("k = %.36f\n", k(a, b, p));
				return 0;
			}
			i++;
			p0 = p;
		}
	}

	return 0;
}
