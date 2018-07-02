Estimativas de máximo para uma série com coeficientes aleatórios de distribuições de máxima entropia.

Compilar o programa que encontra os parâmetros lambda e k:

`gcc -o newton -lm newton.c`

Compilar o programa que gera a função do EP:

`gcc -o lab3 -lm lab3.c`

O programa `newton` recebe quatro parâmetros:

`./newton a b T N`

Onde:

`a`: Início do intervalo

`b`: Fim do intervalo

`T`: Tolerância de aceitação para a raiz

`N`: Número máximo de iterações

O programa `lab3` recebe cinco parâmetros:

`./lab3 a b N M s`

Onde:

`a`: Início do intervalo

`b`: Fim do intervalo

`N`: N pares de coeficientes (fim da somatória)

`M`: Quantidade de valores da função h(x) no grid

`s`: Semente aleatória para gerar os valores uniformes
