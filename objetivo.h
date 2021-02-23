/*
 * objetivo.h
 *
 *  Created on: 03/11/2020
 *
 */

#ifndef TRANFORMADAWD_H_
#define TRANFORMADAWD_H_

#include "mem_structs.h"

void fwht_transform(const int *src, int *dst, int n);
void NoLinealidad(INDIVIDUO *p);
void SAC_0(INDIVIDUO *p);
void aptitud(INDIVIDUO *p);

#endif /* TRANFORMADA_H */
