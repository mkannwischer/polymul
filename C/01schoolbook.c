/**
 * @file 01schoolbook.c
 * @brief Section 2.2.1: Schoolbook multiplication
 *
 * Illustrates schoolbook multiplication in various polynomial rings.
 * This is not rocket science, but for reference of what the other algorithms
 * should be implementing.
 */

#include "poly.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

/**
 * @brief Multiplies two polynomials a and b in Z_q[x]
 *
 * @param c full product with aN+bN-1 coefficients
 * @param a first multiplicand with aN coefficients
 * @param aN number of coefficients in a
 * @param b second multiplicand with bN coefficients
 * @param bN number of coefficients in b
 * @param q modulus
 */
static void polymul_schoolbook(T* c, const T* a, size_t aN, const T* b,
                               size_t bN, T q){
    size_t i,j;
    uint32_t t;
    for(i=0; i<aN+bN-1; i++){
        c[i] = 0;
    }

    for(i=0;i<aN;i++){
        for(j=0;j<bN;j++){
            // since Q is very small, we could compute a 16 bits product only
            t = ((uint32_t)a[i] * b[j]) % q;
            c[i+j] = (c[i+j] + t) % q;
        }
    }
}

/**
 * @brief Cyclic convolution of two polynomials, i.e., multiplication mod x^n - 1
 *
 * @param c product ab mod x^n-1 with n coefficients
 * @param a first multiplicand with n coefficients
 * @param b second multiplicand with n coefficients
 * @param n number of coefficients in a, b, c
 * @param q modulus
 */
static void polymul_schoolbook_cyclic(T* c, const T* a, const T* b,
                               size_t n, T q){
    T t[2*n-1];
    size_t i;
    // perform double-sized schoolbook
    polymul_schoolbook(t, a, n, b, n, q);

    // reduce mod (X^n-1, q)
    for(i=0;i<n;i++){
        c[i] = t[i];
    }
    for(i=0;i<n-1;i++){
        c[i] = (c[i] + t[i+n]) % q;
    }
}

/**
 * @brief Negacyclic convolution of two polynomials, i.e., multiplication mod x^n + 1
 *
 * @param c product ab mod x^n+1 with n coefficients
 * @param a first multiplicand with n coefficients
 * @param b second multiplicand with n coefficients
 * @param n number of coefficients in a, b, c
 * @param q modulus
 */
static void polymul_schoolbook_negacyclic(T* c, const T* a, const T* b,
                               size_t n, T q){
    T t[2*n-1];
    size_t i;
    // perform double-sized schoolbook
    polymul_schoolbook(t, a, n, b, n, q);

    // reduce mod (X^n+1, q)
    for(i=0;i<n;i++){
        c[i] = t[i];
    }
    for(i=0;i<n-1;i++){
        c[i] = (c[i] + q - t[i+n]) % q;
    }
}

int main(void){
    size_t n = 2;
    T q = 17;
    T a[2] = {0, 4};
    T b[2] = {13, 4};

    poly_print("a=", a, n); // a= 4x^1
    poly_print("b=", b, n); // b= 4x^1 + 13x^0

    printf("# Polynomial multiplication in Z_{17}[x]\n");
    T c2n[2*n-1];
    polymul_schoolbook(c2n, a, n, b, n, q);
    poly_print("a*b=", c2n, 2*n-1); // a*b= 16x^2 + 1x^1

    printf("# Polynomial multiplication in Z_{17}[x]/(x^2-1)\n");
    T c[n];
    polymul_schoolbook_cyclic(c, a, b, n, q);
    poly_print("a*b= ", c, n); // a*b= 1x^1 + 16x^0

    printf("# Polynomial multiplication in Z_{17}[x]/(x^2+1)\n");
    polymul_schoolbook_negacyclic(c, a, b, n, q);
    poly_print("a*b= ", c, n ); // a*b= 1x^1 + 1x^0

    return 0;
}
