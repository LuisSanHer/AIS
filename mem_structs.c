/*
 * mem_structs.h
 *
 *  Created on: 08/10/2020
 *
 */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "mem_structs.h"

void alloc_pop(POBLACION *P, int psize){
	size_t i;
	P->size = psize;
	P->ind = (INDIVIDUO*)malloc(sizeof(INDIVIDUO) * psize);
	for(i=0 ; i<psize ; i++){
		P->ind[i].x=(int*)malloc(sizeof(int) * mop.nbin);
		P->ind[i].esp=(int*)malloc(sizeof(int) * mop.nbin);
	}
}

void free_pop(POBLACION *P){
	size_t i;
	size_t psize = P->size;
	for(i=0 ; i< psize ; i++){
		free(P->ind[i].x);
		free(P->ind[i].esp);
	}
	free(P->ind);
}


void cpy_pop(POBLACION *P, POBLACION *Q){
	size_t i;
	assert(P->size == Q->size);
	for(i=0 ; i<P->size ; i++){
		cpy_ind(&P->ind[i], &Q->ind[i]);
	}
}

void cpy_ind(INDIVIDUO *A, INDIVIDUO *B){
	A->f = B->f;
	A->NL = B->NL;
	A->SAC = B->SAC;
	memcpy(A->x, B->x, sizeof(int) * mop.nbin);
	memcpy(A->esp, B->esp, sizeof(int) * mop.nbin);
}
