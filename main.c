#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "rand.h"
#include "ais.h"
#include "objetivo.h"
#include "mem_structs.h"

MOP mop;
AIS ais;
int n;

/*
set title 'Convergencia de la media (n=3)'
set grid
set autoscale
set xlabel 'Generaciones'
set ylabel 'Valor de aptitud'
set key bottom
plot 'MediaN3.txt' u 1:2 w lp lt 1 lw 2 t "Genetico Simple", 'MediaN3.txt' u 1:3 w lp lt 2 lw 3 t "Genetico elitismo", 'MediaN3.txt' u 1:4 w lp lt 3 lw 2 t "Genetico Miu+Lambda", 'MediaN3.txt' u 1:5 w lp lt 4 lw 3 t "SIA"
*/



int main(int argc, char *argv[]){
	srand(time(NULL));

	POBLACION P,Q,Clones, Mejores;
	size_t i=0;
	float semilla;
	int alg, longitud, poblacion, gmax;
	char r[10], p[10], max[10], sem[10];

	printf("Iniciando Sistema inmune artificial...\n");

	do {
	printf("Ingrese exponente (3,5,8,10,12) para la longitud de la cadena binaria (genotipo): ");
	scanf("%s", r);
	longitud = atoi(r);
	//printf("%d\n", longitud);
	if (longitud!=3 && longitud!=5 && longitud!=8 && longitud!=10 && longitud!=12)
	  printf("\tOpción invalida\n");
	} while(longitud!=3 && longitud!=5 && longitud!=8 && longitud!=10 && longitud!=12);


	do {
	  printf("\nIngrese el tamaño de la población: ");
	  scanf("%s", p);
	  poblacion = atoi(p);
	  //printf("%d\n", poblacion);
	  if (poblacion%2 != 0 || poblacion == 0 )
	    printf("\tIngrese un número par\n");
	} while(poblacion%2 != 0 || poblacion == 0);


	do {
	  printf("\nIngrese el máximo de generaciones: ");
	  scanf("%s", max);
	  gmax = atoi(max);
	  //printf("%d\n", gmax);
	  if (gmax == 0 )
	    printf("\tIngrese un número \n");
	} while(gmax == 0);

	do {
	  printf("\nIngrese el una semilla(1-100): ");
	  scanf("%s", sem);
	  semilla = atoi(sem);
	  //printf("%d\n", semilla);
	  if ( semilla<1 || semilla>100 )
	    printf("\tOpción invalida\n");
	} while( semilla<1 || semilla>100 );

	n = longitud;
	mop.nbin = pow(2.0, n);
	//mop.nbin = 500;
	mop.nobj = 1;

	//ais.psize = poblacion;
	ais.psize = poblacion;
	ais.n_mejores = ais.psize*0.30; //30% de la población, para seleccionar mejores
	ais.n_peores = ais.psize*0.40;  //40% de la población, para reemplazar peores
	ais.Gmax = gmax;
	randomize(semilla/100.0);

	alloc_pop(&P, ais.psize);
	alloc_pop(&Mejores, ais.n_mejores);
	FILE* file = fopen("Experimentos/experimento.txt", "w");
	Inicializar(&P); //Generar población inicial de anticuerpos aleatoriamente.
	for(i=0 ; i<ais.Gmax ; i++){
	  Evaluacion(&P);  //Presentar el antígeno y calcular la afinidad del anticuerpo por el antígeno.
	  Seleccionar(&P, &Mejores, ais.n_mejores); //Seleccionar los n mejores anticuerpos.
	  Clonacion(&Mejores, &Clones); //Clonación directamente proporcional a la aptitud.
	  Hipermutacion(&Clones); //Mutación inversamente proporcional a la aptitud.
		Autorregulacion(&P, &Q, &Clones); //Unir P con Clones, ordenar por aptitud y reemplazar los peores por aleatorios
		Estadisticas(&P, i, file);

		free_pop(&Q);
		free_pop(&Clones);
	}
	fclose(file);
	free_pop(&Mejores);
	free_pop(&P);

}
