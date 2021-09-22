/**
 * @file poly.h
 * @brief Common code for working with polynomials.
 *
 * Polynomials are represented as T[] with T being defined in common.h.
 */
#ifndef POLY_H
#define POLY_H
#include "common.h"
#include <stdint.h>
#include <stddef.h>

/**
 * @brief Dumps a polynomial to stdout, e.g., "a= 16x^2 + 1x^1".
 *
 * @param polyname name to be printed, e.g, "a" in the example above
 * @param a polynomial
 * @param n number of coefficients
 */
void poly_print(char *polyname, T *a, size_t n);

/**
 * @brief Dumps a polynomial to stdout, e.g., [0,1,16].
 *
 * @param a polynomial
 * @param n number of coefficients
 */
void poly_print2(T *a, size_t n);

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
void poly_polymul_ref(T *c, const T *a, size_t aN, const T *b, size_t bN, T q);

/**
 * @brief Samples a uniformly random polynomial with coefficients in [0, q).
 *
 * @param a polynomial
 * @param n number of coefficients in polynomial
 * @param q modulus
 */
void poly_random(T *a, size_t n, T q);

/**
 * @brief Compares two polynomials mod Q.
 *
 * @param a first polynomial
 * @param b second polynomial
 * @param n number of coefficients in a and b
 * @param q modulus
 * @return int 0 if coefficients are equal mod q, 1 otherwise
 */
int poly_compare(T *a, T *b, size_t n, T q);

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
void polymul_cyclic_ntt_forward_reference(T *t, const T *a, T q, T n, T root);

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
void polymul_cyclic_ntt_inverse_reference(T *t, const T *a, T q, T n, T root);

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
void polymul_negacyclic_ntt_forward_reference(T *t, const T *a, T q, T n, T root);

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
void polymul_negacyclic_ntt_inverse_reference(T *t, const T *a, T q, T n, T root);
#endif
