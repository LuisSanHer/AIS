/*
 * ais.c
 *
 *  Created on: 20/01/2021
 *
 */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "objetivo.h"
#include "rand.h"
#include "ais.h"

void Inicializar(POBLACION *P){
	size_t i,j;
	for(i=0 ; i< P->size ; i++){
		for(j=0 ; j<mop.nbin ; j++){
			if(randomperc() < 0.5){
				P->ind[i].x[j] = 0;
			}else{
				P->ind[i].x[j] = 1;
			}
		}
	}
}

void Evaluacion(POBLACION *Q){
	size_t i;
	int n = Q->size;
	for(i=0 ; i<n ; i++){ //Para todos los individuos
		NoLinealidad(&Q->ind[i]); //Calcular la NoLinealidad.
		SAC_0(&Q->ind[i]);				//Calcular el criterio de avalancha estricto.
		aptitud(&Q->ind[i]);      //Calcular aptitud respecto al SAC y NL (restricciones).

	}
}

void Seleccionar(POBLACION *P, POBLACION *Mejores, int n){
	Ordenar(P);
	size_t i;
	for(i=0 ; i<n ; i++){
		cpy_ind(&Mejores->ind[i], &P->ind[i]);
	}
}

void Clonacion(POBLACION *Mejores, POBLACION *Clones){
  int count=0, k=0;
	// Calcular cuantos clones serán creados
	int tam_arr_clones = Calcular_clones(Mejores, ais.n_mejores);
	// Reservar memoria para clones.
	alloc_pop(Clones, tam_arr_clones);

  float best_afin = Mejores->ind[0].f, aux_afin=0.0;
  for (size_t i=0 ; i<Mejores->size ; i++) {
    aux_afin = Mejores->ind[i].f;
    count = (ais.psize/2)/(aux_afin-best_afin+1);
    for (size_t j=0 ; j<count ; j++) {
      cpy_ind(&Clones->ind[k], &Mejores->ind[i]);
      k++;
    }
  }
}

void Hipermutacion(POBLACION *Clones){
	size_t i, n = Clones->size;
  double Pm = 0.0;
  double worst_afin = Clones->ind[n-1].f;
	double best_afin = Clones->ind[0].f;
	double aux_afin=0.0;
	for(i=0 ; i<n ; i++){ //Para todos los individuos
		aux_afin = Clones->ind[i].f;
		Pm = (aux_afin*(1.0/mop.nbin))/best_afin;
		bit_wise_mutation(&Clones->ind[i], Pm); //Mutar bit a bit
	}
}

void bit_wise_mutation(INDIVIDUO *A, double Pm){
	size_t i;
	for(i=0 ; i<mop.nbin ; i++){ //Para toda la cadena binaria
		if(randomperc() < Pm){ //Si random caé dentro del Parametro de mutación
			A->x[i] = 1 - A->x[i]; // MUTAR (Cambia de 0 a 1 y viceversa).
		}
	}
}

void Reemplazar(POBLACION *P, int n){
	size_t x = P->size - n;
	for (size_t i=x ; i<P->size ; i++) {
		for(size_t j=0 ; j<mop.nbin ; j++){ //Generar nuevo
			if(randomperc() < 0.5){           //individio
				P->ind[i].x[j] = 0;             //de forma
			}else{                            //aleatoria.
				P->ind[i].x[j] = 1;
			}
		}
	}
}

int Calcular_clones(POBLACION *P, int n){
  int cant_clones=0, tmp = 0;
  float best_afin = P->ind[0].f, aux_afin=0.0;
  for (size_t i=0 ; i<n ; i++) {
    aux_afin = P->ind[i].f;
    tmp = (ais.psize/2)/(aux_afin-best_afin+1);
    cant_clones += tmp;
  }
  return cant_clones;
}

void Display_ind(INDIVIDUO ind){
	size_t i;
	printf("  \033[1;41m NL: %.3lf\033[0m", ind.NL);
	printf("  \033[1;41m SAC: %.3lf\033[0m", ind.SAC);
	printf("  \033[1;41m f: %.3lf\033[0m x: ", - ind.f);
	for(i=0 ; i<mop.nbin ; i++){
		printf("%d", ind.x[i]);
	}
	mop.nbin >= 121?printf("\n"):printf("\n\n");
}

void Display_pop(POBLACION *P){
	size_t i,j;
	for(i=0 ; i<P->size ; i++){
  	printf("Individuo %zu: f: %lf \t x: ", i, P->ind[i].f);
  	for(j=0 ; j<mop.nbin ; j++){
  		printf("%d", P->ind[i].x[j]);
  	}
  	printf("\n");
	}
}

int Mejor_solucion(POBLACION *P){
	size_t i, index;
  INDIVIDUO *mejor = (INDIVIDUO*)malloc(sizeof(INDIVIDUO));
  mejor->x=(int*)malloc(sizeof(int) * mop.nbin);
	mejor->esp=(int*)malloc(sizeof(int) * mop.nbin);
  cpy_ind(mejor, &P->ind[0]);
	index = 0;
	for(i=0 ; i<P->size ; i++){
		if( mejor->f > P->ind[i].f ){	// Evaluacion en términos de
			cpy_ind(mejor, &P->ind[i]);// minimización.
			index = i;                  //
		}
	}
	free(mejor->esp);
  free(mejor->x);
  free(mejor);
	return index;
}

void Estadisticas(POBLACION *P, size_t i, FILE* file){
	int mejor;
	printf("\033[1;44mGeneración: %.3zu\033[0m", i);
	mejor = Mejor_solucion(P);
	Display_ind(P->ind[mejor]);
	fprintf(file,"%lf\n", P->ind[mejor].f);
}

int compare(const void *_a, const void *_b){
	INDIVIDUO *a = (INDIVIDUO *) _a;
	INDIVIDUO *b = (INDIVIDUO *) _b;
	if(a->f < b->f)
		return 0;
	else
		return 1;
}

void Ordenar(POBLACION *T) {
	qsort(T->ind,T->size,sizeof(INDIVIDUO),&compare);
}

void Unir_poblaciones(POBLACION *P, POBLACION *Clones, POBLACION *Q){
	size_t i;
	//Tamaño de la poblacion original + clones
	int tam = P->size + Clones->size;
	//Reservar memoria para unir conjuntos
	alloc_pop(Q, tam);
	for(i=0 ; i<Q->size ; i++){
		if(i < P->size){
			cpy_ind(&Q->ind[i], &P->ind[i]);
		}else{
			cpy_ind(&Q->ind[i], &Clones->ind[i-P->size]);
		}
	}
  Ordenar(Q);
}
