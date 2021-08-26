#include "poly.h"
#include "randombytes.h"
#include <stdio.h>
#include "zq.h"
#include <string.h>
// TODO: test with 32-bit coeffs

void poly_print(char *polyname, T *a, size_t n)
{
  int i;
  int first = 1;
  printf("%s=\t", polyname);
  for (i = n - 1; i >= 0; i--)
  {
    if (a[i] == 0)
      continue;
    if (first)
    {
      first = 0;
    }
    else
    {
      printf(" + ");
    }
    printf("%dx^%d", a[i], i);
  }
  if (first)
    printf("0");
  printf("\n");
}

void poly_print2(T *a, size_t n)
{
  size_t i;
  int first = 1;
  for (i = 0; i < n; i++)
  {
    if (first)
    {
      printf("[");
      first = 0;
    }
    else
    {
      printf(",");
    }
    printf("%d", a[i]);
  }
  printf("]\n");
}

// performs a schoolbook multiplication computing the full product with 2*n-1 coeffs
void poly_polymul_ref(T *c, const T *a, size_t aN, const T *b, size_t bN, 
                             T q)
{
  size_t cN = aN + bN - 1;

  size_t i, j;
  T2 t;
  for (i = 0; i < cN; i++)
  {
    c[i] = 0;
  }

  for (i = 0; i < aN; i++)
  {
    for (j = 0; j < bN; j++)
    {
      t = zq_mod((T2)a[i] * b[j], q);
      c[i + j] = zq_mod(c[i + j] + t, q);
    }
  }
}

// samples a uniformly random polynomial with coefficients in [0, Q)
void poly_random(T *a, size_t n, T q)
{
  size_t i;
  for (i = 0; i < n; i++)
  {
    do
    {
      randombytes((uint8_t *)&a[i], sizeof(T));
    } while (a[i] >= q);
  }
}

// compares two polynomials mod Q
int poly_compare(T *a, T *b, size_t n, T q)
{
  size_t i;
  for (i = 0; i < n; i++)
  {
    if (a[i] % q != b[i] % q)
    {
      return 1;
    }
  }
  return 0;
}

void polymul_cyclic_ntt_forward_reference(T *t, const T *a, T q, T n, T root){
    T i, j;
    T tmp;
    T acopy[n]; // allow inplace transform
    memcpy(acopy, a, sizeof(acopy));

    // init output to zero
    for(i=0;i<n;i++) {
        t[i] = 0;
    }

    // compute the NTT
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            tmp =  ((T2)acopy[j]*zq_pow(root, (uint32_t)i*j, q)) % q;
            t[i] = (t[i] + tmp) % q;
        }
    }
}

void polymul_cyclic_ntt_inverse_reference(T *t, const T *a, T q, T n, T root){
    T i, j;
    T tmp;
    T acopy[n]; // allow inplace transform
    memcpy(acopy, a, sizeof(acopy));
    // init output to zero
    for(i=0;i<n;i++){
        t[i] = 0;
    }

    T rootinv = zq_inverse(root, q);
    // compute the inv ntt
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            tmp = ((T2) acopy[j]*zq_pow(rootinv, (uint32_t)i*j, q))%q;
            t[i] = (t[i] + tmp) % q;
        }
    }

    // multiply by n^-1
    T ninv = zq_pow(n, q-2, q);
    for(i=0;i<n;i++){
        t[i] = ((T2)t[i]*ninv)%q;
    }
}


void polymul_negacyclic_ntt_forward_reference(T *t, const T *a, T q, T n, T root){
    T i, j;
    T tmp;
    T acopy[n]; // allow inplace transform
    memcpy(acopy, a, sizeof(acopy));
    // init output to zero
    for(i=0;i<n;i++) {
        t[i] = 0;
    }

    // compute the NTT
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            tmp = ((T2)acopy[j]*zq_pow(root, (uint32_t)2*i*j +j, q)) % q;
            t[i] = (t[i] + tmp) % q;
        }
    }
}

void polymul_negacyclic_ntt_inverse_reference(T *t, const T *a, T q, T n, T root){
    T i, j;
    T tmp;
    T acopy[n]; // allow inplace transform
    memcpy(acopy, a, sizeof(acopy));
    // init output to zero
    for(i=0;i<n;i++){
        t[i] = 0;
    }

    T rootinv = zq_inverse(root, q);
    // compute the inv ntt
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            tmp = ((T2) acopy[j]*zq_pow(rootinv, (uint32_t) 2*i*j, q))%q;
            t[i] = (t[i] + tmp) % q;
        }
    }

    // multiply by n^-1 and psi^-i
    T ninv = zq_pow(n, q-2, q);
    for(i=0;i<n;i++){
        t[i] = ((T2)t[i]*ninv) % q;
        t[i] = ((T2)t[i]*zq_pow(rootinv, i, q)) % q;
    }
}
