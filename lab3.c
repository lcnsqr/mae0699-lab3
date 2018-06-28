#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Valor uniformemente distribuído em [0,1)
#define RAND ((double)(rand() >> 1)/((RAND_MAX >> 1) + 1))

// Intervalos das distribuições de máxima entropia 
// para os coeficientes A e B do índice 0 ao 9 e 
// suas respectivas constantes lambda e k
struct CoefParam {
	double a, b, l, k;
};
struct CoefParam A[5];
struct CoefParam B[5];

// Parâmetros das exponenciais
double expParam[5];

// Variável aleatória exponencial
double exponencial(double lambda){
	return log(1-RAND)/(-lambda);
}

// Variável aleatória normal (Box-Muller)
double normal(double mu, double sigma){
	double u1, u2, w, m;
	static double x1, x2;
	static int call = 0;

	if (call == 1) {
		// Retornar a variável gerada anteriormente
		call = !call;
		return mu + sigma * x2;
	}

	// Aceitação/rejeição dentro do círculo unitário
	do {
		u1 = -1.0 + RAND * 2.0;
		u2 = -1.0 + RAND * 2.0;
		w = pow (u1, 2) + pow (u2, 2);
	} while (w >= 1 || w == 0);

	// Fator comum
	m = sqrt((-2 * log (w)) / w);
	x1 = u1 * m;
	x2 = u2 * m;

	call = !call;

	return mu + sigma * x1;
}

// Variável aleatória bernoulli
int bernoulli(double p){
	return ( RAND < p ) ? 1 : 0;
}

// Função de máxima entropia com as constante k, l
double f(double l, double k, double x){
	return k*exp(l*x);
}

// Variável aleatória para a função de máxima entropia no intervalo [a,b]
double rand_f(double a, double b, double l, double k){
	// O valor máximo da função está no menor 
	// extremo do intervalor em relação à origem
	double ymax;
	if ( fabs(a) < fabs(b) ){
		ymax = f(l, k, a);
	}
	else {
		ymax = f(l, k, b);
	}
	// Gerar um ponto uniformemente distribuído no semi-plano [a,b]x[0,max]
	double x, y;
	do {
		x = a + (b-a) * RAND;
		y = ymax * RAND;
		// Aceitação/rejeição na distribuição de máxima entropia
	} while (y > f(l, k, x));
	return x;
}

// Coeficientes "A" para a série
double A_k(int k){
	// Valor do coeficiente
	double c;
	if ( k >= 0 && k <= 4 ){
		// ~ Máxima entropia
		c = rand_f(A[k].a, A[k].b, A[k].l, A[k].k);
	}
	else if (k >= 5 && k <= 9 ){
		// ~ Exponencial com (-1)^(~ bernoulli)
		c = pow(-1, bernoulli(0.5))*exponencial(expParam[k-5]);
	}
	else {
		// Normal se k > 9
		c = normal(0, 2.0/pow(k-9, 4));
	}
	return c;
}

// Coeficientes "B" para a série
double B_k(int k){
	// Valor do coeficiente
	double c;
	if ( k >= 0 && k <= 4 ){
		// ~ Máxima entropia
		c = rand_f(B[k].a, B[k].b, B[k].l, B[k].k);
	}
	else if (k >= 5 && k <= 9 ){
		// ~ Exponencial com (-1)^(~ bernoulli)
		c = pow(-1, bernoulli(0.5))*exponencial(expParam[k-5]);
	}
	else {
		// Normal se k > 9
		c = fabs(normal(0, 4.0/pow(k-9, 4)));
	}
	return c;
}

// Função h(x)
// Calcular série em x até o máximo de iterações N
double h(double x, double *AC, double *BC, int N){
	// Soma
	double s = 0;
	// Contagem da soma
	int k = 0;
	s += AC[k] * sin(2*M_PI*k*x) + BC[k] * cos(2*M_PI*k*x);
	k++;
	while ( k < N ){
		s += AC[k] * sin(2*M_PI*k*x) + BC[k] * cos(2*M_PI*k*x);
		k++;
	}
	// Agregar estimativa do erro de truncamento
	s += normal(0, (2*pow(M_PI, 2))/45.0);
	return s;
}

// Derivada da função h(x)
// Aproximação numérica da derivada usando os pontos do grid
// grid: Conjunto dos valores de h(x)
// k: Aproximação da derivada no k-ésimo ponto
// dx: Espaço entre os pontos do grid
// M: Total de pontos no grid
double dh(double *grid, int k, double dx, int M){
	double d;
	if ( k == 0 ){
		// Primeiro ponto
		d = (-3*grid[k] + 4*grid[k+1] - grid[k+2])/(2*dx);
	}
	else if ( k == M ){
		// Último ponto
		d = (grid[k-2] - 4*grid[k-1] + 3*grid[k])/(2*dx);
	}
	else {
		// Ponto intermediário
		d = (grid[k+1] - grid[k-1])/(2*dx);
	}
	return d;
}

int main(int argc, char **argv){
	// Semente aleatória baseada no último parâmetro da linha de comando
	srand((unsigned int)atoi(argv[6]));

	// Carregar parâmetros de cada coeficiente da série
	A[0].a = -1.0;
	A[0].b = 1.0;
	A[0].l = 0.0;
	A[0].k = 0.5;

	B[0].a = -1.0;
	B[0].b = 0.5;
	B[0].l = 1.432750533271374804300535288348328322;
	B[0].k = 0.792297876791121291617514543759170920;

	A[1].a = -0.5;
	A[1].b = (2.0/3.0);
	A[1].l = -0.743864787874943034218233606225112453;
	A[1].k = 0.883955805769716773667710185691248626;

	B[1].a = (-2.0/3.0);
	B[1].b = 0.5;
	B[1].l = 0.743864787874943145240536068740766495;
	B[1].k = 0.883955805769716884690012648206902668;

	A[2].a = -3.0;
	A[2].b = 7.0;
	A[2].l = -0.267210385527338611932890444222721271;
	A[2].k = 0.128768439536812023815670613657857757;

	B[2].a = -5.0;
	B[2].b = 2.0;
	B[2].l = 0.526768681965973151193338708253577352;
	B[2].k = 0.188402460085155792901545623863057699;

	A[3].a = -0.1;
	A[3].b = 0.9;
	A[3].l = -9.995441133814841450089261343237012625;
	A[3].k = 3.678961817100728559637445869157090783;

	B[3].a = -0.8;
	B[3].b = 0.2;
	B[3].l = 4.801007549722517531165522086666896939;
	B[3].k = 1.853136731909351020419762789970263839;

	A[4].a = -0.2;
	A[4].b = 0.3;
	A[4].l = -2.459866400763915272875692608067765832;
	A[4].k = 2.125242399293599593335102326818741858;

	B[4].a = -0.45;
	B[4].b = 0.05;
	B[4].l = 19.990882267629682900178522686474025249;
	B[4].k = 7.357923634201457119274891738314181566;

	// Parâmetros das exponenciais
	expParam[0] = 1.0;
	expParam[1] = 2.0;
	expParam[2] = 0.5;
	expParam[3] = 3.0;
	expParam[4] = 5.0;

	// Gerar amostras das funções de máxima entropia para os gráficos
	/*
	int i = atoi(argv[1]);
	for (double x = B[i].a; x < B[i].b; x += 1e-2){
		printf("%f %f\n", x, f(B[i].l, B[i].k, x));
	}
	*/

	// Traçar função no intervalo indicado na linha de comando
	double a = atof(argv[1]);
	double b = atof(argv[2]);

	// Armazenar N pares de coeficientes 
	int N = atoi(argv[3]);

	// Quantidade de valores da função h(x) no grid
	int M = atoi(argv[4]);

	// Espaço entre os pontos do grid
	double dx = fabs((b-a)/(double)(M-1));

	// Conjuntos dos coeficientes para a função h(x)
	double *AC = (double*)malloc(N*sizeof(double));
	double *BC = (double*)malloc(N*sizeof(double));

	// Grid com os valores gerados pela função h(x)
	double *grid = malloc(M*sizeof(double));
	// Aproximação numérica da derivada em cada ponto do grid
	double *grid_d = malloc(M*sizeof(double));
	// Estimativa dos valores entre os pontos do grid
	double *grid_m = malloc((M-1)*sizeof(double));

	// Variáveis auxiliares para determinar o encontro das derivadas
	double y[2], m[2], x[2], c[2];

	// Armazenar o maior valor estimado
	double max;

	// Gerar a quantidade solicitada de estimativas para o máximo
	for (int j = 0; j < atoi(argv[5]); j++){
		// Zerar máximo
		max = 0;

		// Gerar e armazenar os coeficientes
		for (int k = 0; k < N; k++){
			AC[k] = A_k(k);
			BC[k] = B_k(k);
		}

		// Computar valores de h(x)
		for (int i = 0; i < M; i++){
			grid[i] = h(a + i*dx, AC, BC, N);
			if ( grid[i] > max ) max = grid[i];
		}

		// Computar derivadas de h(x)
		for (int i = 0; i < M; i++){
			grid_d[i] = dh(grid, i, dx, M-1);
		}

		// Estimar o valor entre os pontos do 
		// grid usando as derivadas em cada ponto.
		for (int i = 0; i < M - 1; i++){
			y[0] = grid[i];
			y[1] = grid[i+1];
			m[0] = grid_d[i];
			m[1] = grid_d[i+1];
			x[0] = a + i*dx;
			x[1] = a + (i+1)*dx;
			c[0] = y[0] - m[0] * x[0];
			c[1] = y[1] - m[1] * x[1];
			if (m[1] - m[0] != 0){
				// Valor estimado entre os pontos i e i+1
				grid_m[i] = (c[0]*m[1] - c[1]*m[0])/(m[1] - m[0]);
			}
			else {
				// Ponto entre y[0] e y[1]
				grid_m[i] = (y[0] + y[1])/2;
			}
			if ( grid_m[i] > max ) max = grid_m[i];
		}

		printf("%f\n", max);
	}

	free(AC);
	free(BC);
	free(grid);
	free(grid_d);
	free(grid_m);
	return 0;
}
