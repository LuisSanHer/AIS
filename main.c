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

int main(int argc, char *argv[]){

	POBLACION P,Q,Clones, Mejores;
	size_t i=0;
	float semilla;
	int alg, longitud, poblacion, gmax;
	char r[10], p[10], max[10], sem[10];

	printf("Iniciando Sistema inmune artificial...\n");

	do {
	printf("Ingrese exponente (4,6,8,10,12) para la longitud de la cadena binaria (genotipo): ");
	scanf("%s", r);
	longitud = atoi(r);
	//printf("%d\n", longitud);
	if (longitud!=4 && longitud!=6 && longitud!=8 && longitud!=10 && longitud!=12)
	  printf("\tOpción invalida\n");
	} while(longitud!=4 && longitud!=6 && longitud!=8 && longitud!=10 && longitud!=12);


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
	  printf("\nIngrese el una semilla(1-10): ");
	  scanf("%s", sem);
	  semilla = atoi(sem);
	  //printf("%d\n", semilla);
	  if ( semilla<1 || semilla>10 )
	    printf("\tOpción invalida\n");
	} while( semilla<0 || semilla>10 );

	n = longitud;
	mop.nbin = pow(2, n);
	//mop.nbin = 500;
	mop.nobj = 1;

	//ais.psize = poblacion;
	ais.psize = poblacion;
	ais.n_mejores = ais.psize*0.25; //25% de la población, para seleccionar mejores
	ais.n_peores = ais.psize*0.25;  //25% de la población, para reemplazar peores
	ais.Gmax = gmax;
	randomize(semilla/10.0);

	alloc_pop(&P, ais.psize);
	alloc_pop(&Mejores, ais.n_mejores);
	FILE* file = fopen("experimento.txt", "w");
	Inicializar(&P); //Generar población inicial de anticuerpos aleatoriamente.
	for(i=0 ; i<ais.Gmax ; i++){
	  Evaluacion(&P);  //Presentar el antígeno y calcular la afinidad del anticuerpo por el antígeno.
	  Seleccionar(&P, &Mejores, ais.n_mejores); //Seleccionar los n mejores anticuerpos.
	  Clonacion(&Mejores, &Clones);
	  //Display_pop(&Clones);
	  Hipermutacion(&Clones);
	  Evaluacion(&Clones);  //Calcular la afinidad vs el antígeno de los anticuerpos clonados y mutados.
	  Unir_poblaciones(&P, &Clones, &Q);
	  Seleccionar(&Q, &P, ais.psize); //Seleccionar los n mejores anticuerpos de Q y guardarlos en P.
		//Display_pop(&Q);
		Estadisticas(&P, i, file);
		Reemplazar(&P, ais.n_peores);
	  //Display_pop(&Q);
		free_pop(&Q);
		free_pop(&Clones);
	}
	fclose(file);
	free_pop(&Mejores);
	free_pop(&P);

}
