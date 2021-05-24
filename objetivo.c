/*  Implementación de la transformada rápida de WH
 * objetivo.c
 *
 *  Created on: 03/11/2020
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "objetivo.h"

int iabs(int a){
  if(a<0)
    return -a;
  return a;
}

double max(double num1, double num2) {
   double result;
   if (num1 > num2)
      result = num1;
   else
      result = num2;
   return result;
}

void fwht_transform(const int *src, int *dst, int n) {
    int *a = (int*)malloc(sizeof(int)*n);
    int *b = (int*)malloc(sizeof(int)*n);
    void *tmp;
    memcpy(a, src, sizeof(int)*n);
    int max_unex=-1;
    // Fast Walsh Hadamard Transform.
    int i, j, s;
    int t;
    for (i = n>>1; i > 0; i>>=1) // i=4  ----->  i<<1 :     00000000000000001000
    {
        for (j = 0; j < n; j++)
        {
            s = j/i%2;
            b[j]=a[(s?-i:0)+j]+(s?-1:1)*a[(s?0:i)+j];
        }
        tmp = a;
        a = b;
        b = tmp;
    }
    memcpy(dst, a, sizeof(int)*n);
    free(a);
    free(b);
}

void NoLinealidad(INDIVIDUO *p){
  int h = 1;
  int x, y;
  fwht_transform(p->x, p->esp, mop.nbin);
  //Obtener el maximo valor absoluto, sin tomar en cuenta el primero
  int mayor = 0;
  for (size_t i=1 ; i<mop.nbin ; i++) {
    if(mayor < iabs(p->esp[i])){
      mayor = iabs(p->esp[i]);
    }
  }
  //CLARK 2005: ECUACION 6
  p->NL = (pow(2, n-1) - mayor); //No linealidad completa sin dividir max R= 112
                               //(CORRECTA) hasta cierto punto...
                               //Es la que se acerca a los resultados de Burnett.
}

void SAC_0(INDIVIDUO *p){
  int w, c, balance;
  //balance
  balance=0;
  for(w=0; w<mop.nbin; w++){
    balance += p->x[w];
  }
  balance = (balance==(mop.nbin/2))? 1:0;
  if(balance==1){
    p->esp[0] = 0;
  }
  else{
    p->esp[0] = p->esp[0]-(mop.nbin/2);
  }
  //aux^2 ;
  int* aux_2 =(int*)malloc(sizeof(int)*mop.nbin);
  for(w=0; w<mop.nbin; w++){
    aux_2[w] = p->esp[w]*p->esp[w];
  }

  int i;
  c = 1; // 00000000001
  int Total, s;
  Total=0;
  for(i=0; i<n; i++){
    s=0;
    for(w=0; w<mop.nbin; w++){
      s += ((w & c)==0)? -aux_2[w] : aux_2[w];
    }
    Total += s;
    c = c<<1;
  }
  //p->SAC = max(0, iabs(Total));
  p->SAC = Total;
  free(aux_2);
}

void aptitud(INDIVIDUO *p){
  double e = 0.0001;
  double v = max(0, (iabs(p->SAC) - e));
  p->f = (-p->NL) + pow(v,2);
}
