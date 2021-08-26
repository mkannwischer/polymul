#ifndef POLY_H
#define POLY_H
#include "common.h"
#include <stdint.h>
#include <stddef.h>


void poly_print(char *polyname, T *a, size_t n);

void poly_print2(T *a, size_t n);

// performs a schoolbook multiplication computing the full product with 2*n-1 coeffs
void poly_polymul_ref(T *c, const T *a, size_t aN, const T *b, size_t bN, T q);

// samples a uniformly random polynomial with coefficients in [0, Q)
void poly_random(T *a, size_t n, T q);

// compares two polynomials mod Q
int poly_compare(T *a, T *b, size_t n, T q);

// just for reference
void polymul_cyclic_ntt_forward_reference(T *t, const T *a, T q, T n, T root);

void polymul_cyclic_ntt_inverse_reference(T *t, const T *a, T q, T n, T root);

void polymul_negacyclic_ntt_forward_reference(T *t, const T *a, T q, T n, T root);

void polymul_negacyclic_ntt_inverse_reference(T *t, const T *a, T q, T n, T root);
#endif