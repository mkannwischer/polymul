/**
 * @file poly.c
 * @brief Common code for working with polynomials.
 *
 * Polynomials are represented as T[] with T being defined in common.h.
 */

#include "poly.h"
#include "randombytes.h"
#include <stdio.h>
#include "zq.h"
#include <string.h>
// TODO: test with 32-bit coeffs

/**
 * @brief Dumps a polynomial to stdout, e.g., "a= 16x^2 + 1x^1".
 *
 * @param polyname name to be printed, e.g, "a" in the example above
 * @param a polynomial
 * @param n number of coefficients
 */
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

/**
 * @brief Dumps a polynomial to stdout, e.g., [0,1,16].
 *
 * @param a polynomial
 * @param n number of coefficients
 */
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

/**
 * @brief Performs a schoolbook multiplication computing the full product with 2*n-1 coeffs.
 *
 * @param c result polynomial
 * @param a first multiplicand polynomial
 * @param aN number of coefficients in a
 * @param b  second multiplicand polynomial
 * @param bN number of coefficients in b
 * @param q modulus
 */
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

/**
 * @brief Samples a uniformly random polynomial with coefficients in [0, q).
 *
 * @param a polynomial
 * @param n number of coefficients in polynomial
 * @param q modulus
 */
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

/**
 * @brief Compares two polynomials mod Q.
 *
 * @param a first polynomial
 * @param b second polynomial
 * @param n number of coefficients in a and b
 * @param q modulus
 * @return int 0 if coefficients are equal mod q, 1 otherwise
 */
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

/**
 * @brief Naive implementation of a cyclic NTT.
 *
 *  Needs an n-th root of unity.
 *  Computes antt_i = sum_j=0^n (a_j  root^(ij)) for 0 <= i < n.
 *
 * @param t output polynomial in NTT domain
 * @param a input polynomial in normal domain
 * @param q modulus
 * @param n number of coefficients in a (and t)
 * @param root n-th root of unity modulo q
 */
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

/**
 * @brief Naive implementation of an inverse cyclic NTT.
 *
 * Needs an n-th root of unity.
 * Computes a_i = 1/n * sum_j=0^n (a_j  root^(-ij)) for 0 <= i < n.
 *
 * @param t output polynomial in normal domain
 * @param a input polynomial in NTT domain
 * @param q modulus
 * @param n number of coefficients in a (and t)
 * @param root n-th root of unity
 */
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

/**
 * @brief Naive implementation of a negacyclic NTT.
 *
 * Needs a 2n-th root of unity.
 * Computes antt_i = sum_j=0^n (a_j  root^(2ij + j)) for 0 <= i < n.
 *
 * @param t output polynomial in NTT domain
 * @param a input polynomial in normal domain
 * @param q modulus
 * @param n number of coefficients in a (and t)
 * @param root 2n-th root of unity
 */
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

/**
 * @brief Naive implementation of an inverse negacyclic NTT.
 *
 * Needs a 2n-th root of unity.
 * Computes a_i = 1/n * root^(-i) * sum_j=0^n (a_j  root^(-2ij)) for 0 <= i < n.
 *
 * @param t output polynomial in normal domain
 * @param a input polynomial in NTT domain
 * @param q modulus
 * @param n number of coefficients in a (and t)
 * @param root 2n-th root of unity
 */
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
