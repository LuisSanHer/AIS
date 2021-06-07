/*
 * ais.h
 *
 *  Created on: 20/01/2021
 *
 */

#ifndef AIS_H_
#define AIS_H_
#include "objetivo.h"
#include "mem_structs.h"

void Inicializar(POBLACION *P);
void Evaluacion(POBLACION *Q);
void Seleccionar(POBLACION *P, POBLACION *Mejores, int n);
void Clonacion(POBLACION *Mejores, POBLACION *Clones);
void Hipermutacion(POBLACION *Clones);
void bit_wise_mutation(INDIVIDUO *A, double Pm);
void Reemplazar(POBLACION *P, int n);
void Autorregulacion(POBLACION *P, POBLACION *Q, POBLACION *Clones);

int Mejor_solucion(POBLACION *P);
void Ordenar(POBLACION *T);
int Calcular_clones(POBLACION *P, int n);
void Unir_poblaciones(POBLACION *P, POBLACION *Clones, POBLACION *Q);

void Estadisticas(POBLACION *P, size_t i, FILE* file);
double rule_3(double aux_afin, double best_afin, char metodo);
void Display_pop(POBLACION *P);
void Display_ind(INDIVIDUO ind);

#endif /* AIS_H */
