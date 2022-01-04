/**
 * @file 04ntt.c
 * @brief Section 2.2.4: Number-theoretic transform
 *
 * Illustrates the cyclic and negacyclic number-theoretic transform.
 * This is a reference implementation that requires time O(n^2) for the transform.
 * For fast implementations, see `05fft.c`
 */
#include "poly.h"
#include "zq.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>
/**
 * @brief Computes forward (2n-1) NTT of n-coefficient polynomial
 *
 * Computes pntt_i = sum_{j=0}^{n} p[j]*root^{ij} for 0 <= i < 2n-1
 *
 * @param t output polynomial (2n-1 coefficients). corresponds to a evaluated at root^0, ..., root^(2n-2), i.e., p transformed into NTT domain
 * @param a input polynomial (n coefficients)
 * @param q modulus
 * @param n number of coefficients in a
 * @param root (2n-1)-th root of unity modulo q
 */
static void polymul_ntt_forward(T *t, const T *a, T q, T n, T root){
    T i, j;
    T tmp;
    // init output to zero
    for(i=0;i<2*n-1;i++) {
        t[i] = 0;
    }

    // compute the NTT
    for(i=0;i<2*n-1;i++){
        for(j=0;j<n;j++){
            tmp =  ((T2)a[j]*zq_pow(root, (uint32_t)i*j, q)) % q;
            t[i] = (t[i] + tmp) % q;
        }
    }
}

/**
 * @brief Computes inverse NTT.
 *
 * Note that the inputs here are in NTT domain, i.e., have (2n-1) coefficients.
 * Computes p_i =  1/n * sum_{j=0}^{2n-1} pntt[j]*root^{-ij} for 0 <= i < 2n-1
 *
 * @param t output polynomial (2n-1 coefficients). a transformed to normal domain
 * @param a input polynomial (2n-1 coefficients)
 * @param q modulus
 * @param n original n; a has 2n-1 coefficients
 * @param root (2n-1)-th root of unity modulo q
 */
static void polymul_ntt_inverse(T *t, const T *a, T q, T n, T root){
    T i, j;
    T tmp;
    // init output to zero
    for(i=0;i<2*n-1;i++){
        t[i] = 0;
    }

    T rootinv = zq_inverse(root, q);
    // compute the inv ntt
    for(i=0;i<2*n-1;i++){
        for(j=0;j<2*n-1;j++){
            tmp = ((T2) a[j]*zq_pow(rootinv, (uint32_t)i*j, q))%q;
            t[i] = (t[i] + tmp) % q;
        }
    }

    // multiply by n^-1
    T ninv = zq_inverse(2*n-1, q);
    for(i=0;i<2*n-1;i++){
        t[i] = ((T2)t[i]*ninv)%q;
    }

}

/**
 * @brief Performs NTT-based multiplication of two polynomials in the ring Zq[x].
 *
 * For computing the product we evaluate the n-coefficient polynomials and b
 * at (2n-1) powers of a (2n-1)-th root of unity, then multiply coefficient wise,
 * and interpolate the product using the inverse NTT.
 * Primitive (2n-1)-th root of unity modulo q needs to exist.
 *
 * @param c output polynomial (2n-1 coefficients)
 * @param a first multiplicand polynomial (n coefficients)
 * @param b second multiplicand polynomial (n coefficients)
 * @param n number of coefficients in a and b
 * @param q modulus
 * @return int 1 if there is an error, 0 otherwise
 */
static int polymul_ntt(T *c, const T *a, const T *b, size_t n, T q){
    if(!zq_isPrime(q)){
        printf("ERROR: ntt requires prime q.\n");
        return 1;
    }

    if(((q-1) % (2*n-1)) != 0){
        printf("ERROR: ntt requires q-1 divisible by (2*n-1).\n");
        return 1;
    }

    // find an (2n-1)th root of unity mod q
    T root = zq_primitiveRootOfUnity(2*n-1, q);
    T antt[2*n-1]; T bntt[2*n-1];

    // compute NTT(a) and NTT(b)
    polymul_ntt_forward(antt, a, q, n, root);
    polymul_ntt_forward(bntt, b, q, n, root);

    // pointwise-multiplication
    // cntt = antt * bntt
    T cntt[2*n-1];
    for(size_t i=0;i<2*n-1;i++){
        cntt[i] = ((T2)antt[i]*bntt[i]) % q;
    }

    // c = NTT^-1(cntt)
    polymul_ntt_inverse(c, cntt, q, n, root);
    return 0;
}

/**
 * @brief Computes forward cyclic n-NTT of n-coefficient polynomial
 *
 * Computes pntt_i = sum_{j=0}^{n} p[j]*root^{ij} for 0 <= i < n
 *
 * @param t output polynomial (n coefficients). corresponds to a evaluated at root^0, ..., root^(n-1), i.e., p transformed into NTT domain
 * @param a input polynomial (n coefficients)
 * @param q modulus
 * @param n number of coefficients in a
 * @param root n-th root of unity modulo q
 */
static void polymul_cyclic_ntt_forward(T *t, const T *a, T q, T n, T root){
    T i, j;
    T tmp;
    // init output to zero
    for(i=0;i<n;i++) {
        t[i] = 0;
    }

    // compute the NTT
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            tmp =  ((T2)a[j]*zq_pow(root, (uint32_t)i*j, q)) % q;
            t[i] = (t[i] + tmp) % q;
        }
    }
}

/**
 * @brief Computes cyclic inverse NTT.
 *
 * Computes p_i =  1/n * sum_{j=0}^{n} pntt[j]*root^{-ij} for 0 <= i < n
 *
 * @param t output polynomial (n coefficients). a transformed to normal domain
 * @param a input polynomial (n coefficients)
 * @param q modulus
 * @param n number of coefficients in a
 * @param root n-th root of unity modulo q
 */
static void polymul_cyclic_ntt_inverse(T *t, const T *a, T q, T n, T root){
    T i, j;
    T tmp;
    // init output to zero
    for(i=0;i<n;i++){
        t[i] = 0;
    }

    // compute the inv ntt
    T rootinv = zq_inverse(root, q);
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            tmp = ((T2) a[j]*zq_pow(rootinv, (uint32_t)i*j, q))%q;
            t[i] = (t[i] + tmp) % q;
        }
    }

    // multiply by n^-1
    T ninv = zq_inverse(n, q);
    for(i=0;i<n;i++){
        t[i] = ((T2)t[i]*ninv)%q;
    }

}

/**
 * @brief Performs NTT-based multiplication of two polynomials in the ring Zq[x]/(x^n-1), i.e., cyclic.
 *
 * For computing the product we evaluate the n-coefficient polynomials and b
 * at n powers of a n-th root of unity, then multiply coefficient wise,
 * and interpolate the product using the inverse NTT.
 * Primitive n-th root of unity modulo q needs to exist.
 *
 * @param c output polynomial (n coefficients)
 * @param a first multiplicand polynomial (n coefficients)
 * @param b second multiplicand polynomial (n coefficients)
 * @param n number of coefficients in a and b
 * @param q modulus
 * @return int 1 if there is an error, 0 otherwise
 */
static int polymul_cyclic_ntt(T *c, const T *a, const T *b, size_t n, T q){
    if(!zq_isPrime(q)){
        printf("ERROR: cyclic ntt requires prime q.\n");
        return 1;
    }

    if(((q-1) % n) != 0){
        printf("ERROR: cyclic ntt requires q-1 divisible by n\n");
        return 1;
    }

    // find an n-th root of unity mod q
    T root = zq_primitiveRootOfUnity(n, q);
    T antt[n]; T bntt[n];

    // compute NTT(a) and NTT(b)
    polymul_cyclic_ntt_forward(antt, a, q, n, root);
    polymul_cyclic_ntt_forward(bntt, b, q, n, root);

    // pointwise-multiplication
    // cntt = antt * bntt
    T cntt[n];
    for(size_t i=0;i<n;i++){
        cntt[i] = ((uint32_t)antt[i]*bntt[i]) % q;
    }

    // c = NTT^-1(cntt)
    polymul_cyclic_ntt_inverse(c, cntt, q, n, root);
    return 0;
}

/**
 * @brief Compute forward negacyclic n-NTT of n-coefficient polynomial.
 *
 * Computes pntt_i = sum_{j=0}^{n} p[j]*root^{j}*root^{2ij} for 0 <= i < n.
 * This is the same as first twisting to Zq[x]/(x^n-1) and then performing a cyclic NTT.
 *
 * @param t output polynomial (n coefficients)
 * @param a input polynomial (n coefficients)
 * @param q modulus
 * @param n number of coefficients in a
 * @param root 2n-th root of unity modulo q
 */
static void polymul_negacyclic_ntt_forward(T *t, const T *a, T q, T n, T root){
    T i, j;
    T tmp;
    // init output to zero
    for(i=0;i<n;i++) {
        t[i] = 0;
    }

    // compute the NTT
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            tmp = ((T2)a[j]*zq_pow(root, (uint32_t) 2*i*j+j, q)) % q;
            t[i] = (t[i] + tmp) % q;
        }
    }
}

/**
 * @brief Computes negacyclic inverse NTT.
 *
 * Compute p_i =  1/n * root^{-i} * sum_{j=0}^{n} pntt[j]*root^{-2ij} for 0 <= i < n
 *
 * @param t output polynomial (n coefficients)
 * @param a input polynomial (n coefficients)
 * @param q modulus
 * @param n number of coefficients in a
 * @param root 2n-th root of unity modulo q
 */
static void polymul_negacyclic_ntt_inverse(T *t, const T *a, T q, T n, T root){
    T i, j;
    T tmp;
    // init output to zero
    for(i=0;i<n;i++){
        t[i] = 0;
    }

    T rootinv = zq_inverse(root, q);
    // compute the inv ntt
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            tmp = ((T2) a[j]*zq_pow(rootinv, (uint32_t) 2*i*j, q))%q;
            t[i] = (t[i] + tmp) % q;
        }
    }

    // multiply by n^-1 and root^-i
    T ninv = zq_inverse(n, q);
    for(i=0;i<n;i++){
        t[i] = ((T2)t[i]*ninv) % q;
        t[i] = ((T2)t[i]*zq_pow(rootinv, i, q)) % q;
    }
}

/**
 * @brief Performs NTT-based multiplication of two polynomials in the ring Zq[x]/(x^n+1), i.e., negacyclic.
 *
 * For computing the product, we first twist a and b to Zq[x]/(x^n-1) by
 * multiplying by powers of the 2n-th root of unity,
 * then evaluate at n powers of a n-th root of unity,
 * then multiply coefficient-wise,
 * then interpolate the product using the inverse NTT,
 * then twist back to Zq[x]/(x^n+1) by multiplying by power of the inverse
 * 2n-th root of unity.
 * Primitive 2n-th (and consequently, n-th) root of unity modulo q needs to exist.
 *
 * @param c output polynomial (n coefficients)
 * @param a first multiplicand polynomial (n coefficients)
 * @param b second multiplicand polynomial (n coefficients)
 * @param n number of coefficients in a and b
 * @param q modulus
 * @return int 1 if there is an error, 0 otherwise
 */
static int polymul_negacyclic_ntt(T *c, const T *a, const T *b, size_t n, T q){
    if(!zq_isPrime(q)){
        printf("ERROR: negacyclic ntt requires prime q.\n");
        return 1;
    }

    if(((q-1) % (2*n)) != 0){
        printf("ERROR: negacyclic ntt requires q-1 divisible by 2*n\n");
        return 1;
    }
    // find an 2n-th root of unity mod q
    T root = zq_primitiveRootOfUnity(2*n, q);
    T antt[n]; T bntt[n];

    // compute NTT(a) and NTT(b)
    polymul_negacyclic_ntt_forward(antt, a, q, n, root);
    polymul_negacyclic_ntt_forward(bntt, b, q, n, root);

    // pointwise-multiplication
    // cntt = antt * bntt
    T cntt[n];
    for(size_t i=0;i<n;i++){
        cntt[i] = ((uint32_t)antt[i]*bntt[i]) % q;
    }

    // c = NTT^-1(cntt)
    polymul_negacyclic_ntt_inverse(c, cntt, q, n, root);
    return 0;
}

/**
 * @brief Random test of `polymul_ntt`
 *
 * @param n number of coefficients of input polynomials
 * @param q modulus
 * @param printPoly flag for printing inputs and outputs
 * @return int  0 if test is successful, 1 otherwise
 */
static int testcase_ntt(size_t n, T q, int printPoly){
    int rc = 0;
    T a[n], b[n];
    T c_ref[2*n-1], c[2*n-1];

    printf("Testing naive NTT multiplication with q=%d, n=%zu\n", q, n);

    poly_random(a, n, q);
    poly_random(b, n, q);
    if(printPoly) poly_print("a", a, n);
    if(printPoly) poly_print("b", b, n);
    // compute reference product for comparison
    poly_polymul_ref(c_ref, a, n, b, n, q);
    if(printPoly) poly_print("c_ref", c_ref, 2*n-1);

    rc |= polymul_ntt(c, a, b, n, q);
    if(printPoly) poly_print("c", c, 2*n-1);

    if(poly_compare(c_ref, c, 2*n-1, q)) rc = 1;
    return rc;
}

/**
 * @brief Random test of `polymul_cyclic_ntt`
 *
 * @param n number of coefficients of input polynomials
 * @param q modulus
 * @param printPoly flag for printing inputs and outputs
 * @return int  0 if test is successful, 1 otherwise
 */
static int testcase_cyclic_ntt(size_t n, T q, int printPoly){
    int rc = 0;
    T a[n], b[n];
    T c_ref[2*n-1], c[n];

    printf("Testing naive cyclic NTT multiplication with q=%d, n=%zu\n", q, n);

    poly_random(a, n, q);
    poly_random(b, n, q);
    if(printPoly) poly_print("a", a, n);
    if(printPoly) poly_print("b", b, n);
    // compute reference product for comparison
    poly_polymul_ref(c_ref, a, n, b, n, q);
    // reduce reference result mod X^n - 1
    for(size_t i = 0; i < n-1; i++){
        c_ref[i] = (c_ref[i] + c_ref[n+i]) % q;
    }
    if(printPoly) poly_print("c_ref", c_ref, n);

    rc |= polymul_cyclic_ntt(c, a, b, n, q);
    if(printPoly) poly_print("c", c, n);

    if(poly_compare(c_ref, c, n, q)) rc = 1;
    return rc;
}

/**
 * @brief Random test of `polymul_negacyclic_ntt`
 *
 * @param n number of coefficients of input polynomials
 * @param q modulus
 * @param printPoly flag for printing inputs and outputs
 * @return int  0 if test is successful, 1 otherwise
 */
static int testcase_negacyclic_ntt(size_t n, T q, int printPoly){
    int rc = 0;
    T a[n], b[n];
    T c_ref[2*n-1], c[n];

    printf("Testing naive negacyclic NTT multiplication with q=%d, n=%zu\n", q, n);

    poly_random(a, n, q);
    poly_random(b, n, q);
    if(printPoly) poly_print("a", a, n);
    if(printPoly) poly_print("b", b, n);
    // compute reference product for comparison
    poly_polymul_ref(c_ref, a, n, b, n, q);
    // reduce reference result mod X^n - 1
    for(size_t i = 0; i < n-1; i++){
        c_ref[i] = (c_ref[i] + q - c_ref[n+i]) % q;
    }
    if(printPoly) poly_print("c_ref", c_ref, n);

    rc |= polymul_negacyclic_ntt(c, a, b, n, q);
    if(printPoly) poly_print("c", c, n);

    if(poly_compare(c_ref, c, n, q)) rc = 1;
    return rc;
}

int main(void){
    int rc = 0;

    // plain NTT (no convolution)
    rc |= testcase_ntt(4, 29, 1);
    rc |= testcase_ntt(256, 3067, 0);

    // cyclic NTT (mod x^n-1)
    rc |= testcase_cyclic_ntt(8, 17, 1);
    rc |= testcase_cyclic_ntt(256, 3329, 0);

    // negacyclic NTT (mod x^n+1)
    rc |= testcase_negacyclic_ntt(8, 17, 1);
    rc |= testcase_negacyclic_ntt(256, 7681, 0);

    if(rc){
        printf("ERROR.\n");
    } else {
        printf("ALL GOOD.\n");
    }

    return rc;
}
