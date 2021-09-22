/**
 * @file 03toom.c
 * @brief Section 2.2.3: Toom--Cook multiplication
 *
 * Illustrates Toom-3, and Toom-4 multiplication for any modulus.
 */
#include "poly.h"
#include "zq.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

/**
 * @brief Compute polynomial product using three-way Toom--Cook multiplication
 * Assumes n is the same for both a and b; also 3 must divide n.
 * Sets y=x^(n/3), s.t.
 * a = a0 + y a1 + y^2 a2
 * b = b0 + y b1 + y^2 b2
 *
 * Then evaluates at y = {0, infty, 1, -1, -2}, performs 5 multiplications
 * of degree n/3 polynomials, and interpolates c using the formulas presented
 * in the thesis.
 *
 * During the interpolation, we need to divide by 2 and 3
 * In case q is co-prime to 2 and 3, we can simply multiply by the inverse
 *
 * As q is often a power of two, the inverse of two does often not exist.
 * In that case there is a workaround: We do the small multiplications
 * modulo 2q instead of q, such that we can divide by two during the
 * interpolation without changing the result.
 *
 * @param c product ab with 2n-1 coefficients
 * @param a first multiplicand with n coefficients
 * @param b second multiplicand with n coefficients
 * @param n number of coefficients in a and b. needs to be divisble by 3
 * @param q modulus. needs to be co-prime to 3
 */
static void polymul_toom3(T* c, const T* a, const T* b, size_t n, T q){
    if(n%3 != 0) {
        printf("ERROR: toom3 currently only supports n divisible by 3. please pad your polynomial.");
        return;
    }
    if(q%3 == 0) {
        printf("ERROR: toom3 currently only supports q co-prime to 3.");
        return;
    }
    if(q>(1<<(sizeof(T)*8-1))){
        printf("ERROR: toom3 currently requires q to be at least one bit smaller than the word size.");
        return;
    }

    // During the interpolation, we need to divide by 2 and 3
    // In case q is co-prime to 2 and 3, we can simply multiply by the inverse

    // As q is often a power of two, the inverse of two does often not exist.
    // In that case there is a workaround: We do the small multiplications
    // modulu 2q instead of q, such that we can divide by two during the
    // interpolation without changing the result.

    T inv2 = zq_inverse(2, q);
    T2 qq = q;

    if(inv2 == 0){
        qq = 2*q;
    }

    T inv3 = zq_inverse(3, qq);
    size_t nsplit = n/3;


    // split a and b in 3 parts each, i.e.,
    // a = a0 + y a1 + y^2 a2
    // b = b0 + y b1 + y^2 b2
    const T *a0 = &a[0*nsplit];
    const T *a1 = &a[1*nsplit];
    const T *a2 = &a[2*nsplit];

    const T* b0 = &b[0*nsplit];
    const T *b1 = &b[1*nsplit];
    const T *b2 = &b[2*nsplit];

    // Evaluate at y={0, \infty, 1, -1, -2}
    const T *a_0 = a0; const T *a_inf = a2;
    T a_1[nsplit], a_m1[nsplit], a_m2[nsplit];
    const T *b_0 = b0; const T *b_inf = b2;
    T b_1[nsplit], b_m1[nsplit], b_m2[nsplit];

    for(size_t i=0; i<nsplit; i++){
        //a_1   = a0 + a1 + a2
        //a_m1  = a0 - a1 + a2
        a_1[i] = zq_mod(a0[i] + a2[i], qq);
        a_m1[i]= zq_mod(a_1[i] + qq - a1[i], qq);
        a_1[i] = zq_mod(a_1[i] + a1[i], qq);

        //a_m2  = a0 - 2*a1 + 4*a2
        a_m2[i] = zq_mod((a0[i] + 2*qq - 2*a1[i] + 4*a2[i]), qq);

        //b_1   = b0 + b1 + b2
        //b_m1  = b0 - b1 + b2
        b_1[i] = zq_mod(b0[i] + b2[i], qq);
        b_m1[i]= zq_mod(b_1[i] + qq - b1[i], qq);
        b_1[i] = zq_mod(b_1[i] + b1[i], qq);

        //b_m2  = b0 - 2*b1 + 4*b2
        b_m2[i] = zq_mod(b0[i] + 2*qq - 2*b1[i] + 4*b2[i], qq);
    }


    // do smaller multiplications (in this case using schoolbook; could also be
    // other Toom or Karatsuba)
    T c_0[2*nsplit-1], c_1[2*nsplit-1], c_m1[2*nsplit-1], c_m2[2*nsplit-1],
      c_inf[2*nsplit-1];

    poly_polymul_ref(c_0, a_0, nsplit, b_0, nsplit, qq);
    poly_polymul_ref(c_inf, a_inf, nsplit, b_inf, nsplit, qq);
    poly_polymul_ref(c_1,  a_1, nsplit,  b_1, nsplit, qq);
    poly_polymul_ref(c_m1, a_m1, nsplit, b_m1, nsplit, qq);
    poly_polymul_ref(c_m2, a_m2, nsplit, b_m2, nsplit, qq);


    // interpolation of c = c_0 + v1 y^1 + v2 y^2 + v3 y^3 + c_inf y^4
    for(size_t i=0; i<2*n-1; i++){
        c[i] = 0;
    }

    T v1, v2, v3;
    for(size_t i=0; i<2*nsplit-1; i++){
        // -v1 + v2 - 3v3 + 5c_inf
        v3 = zq_mod((c_m2[i] + qq - c_1[i])*inv3, qq);
        // v1 + v3
        if(inv2){
            v1 = zq_mod((c_1[i] + qq - c_m1[i])*inv2, qq);
        } else {
            v1 = zq_mod((c_1[i] + qq - c_m1[i])>>1, qq);
        }
        // -v1 + v2 -  v3 +  c_inf
        v2 = zq_mod(c_m1[i] + qq - c_0[i], qq);

        // v3
        if(inv2){
            // if 2 is invertible, we multiply by the inverse
            v3 = zq_mod((v2 + qq - v3)*inv2 + 2*c_inf[i], qq);
        } else {
            // otherwise, we divide by two (in which case the small
            // multiplications need to be correct mod 2*q)
            v3 = zq_mod(((v2 + qq - v3)>>1) + 2*c_inf[i], qq);
        }
        // v2
        v2 = zq_mod(v2 + v1  + qq - c_inf[i], qq);
        // v1
        v1 = zq_mod(v1 + qq - v3, qq);

        c[0*nsplit + i] = zq_mod(c[0*nsplit + i]+c_0[i], q);
        c[1*nsplit + i] = zq_mod(c[1*nsplit + i]+v1, q);
        c[2*nsplit + i] = zq_mod(c[2*nsplit + i]+v2, q);
        c[3*nsplit + i] = zq_mod(c[3*nsplit + i]+v3, q);
        c[4*nsplit + i] = zq_mod(c[4*nsplit + i]+c_inf[i], q);
    }
}

/**
 * @brief Compute polynomial product using four-way Toom--Cook multiplication.
 * Assumes n is the same for both a and b; also 4 must divide n.
 * Sets y=x^(n/4), s.t.
 * a = a0 + y a1 + y^2 a2 + y^3 a3
 * b = b0 + y b1 + y^2 b2 + y^3 a3
 *
 * Then evaluates at y = {0, infty, 1, -1, 2, -2, 3}, performs 7 multiplications
 * of degree n/4 polynomials, and interpolates c using the formulas presented
 * in the thesis.
 *
 * During the interpolation, we need to divide by 8 and 3
 * In case q is co-prime to 8 and 3, we can simply multiply by the inverse.
 * As q is often a power of two, the inverse of 8 does often not exist.
 * In that case there is a workaround: We do the small multiplications
 * modulo 8q instead of q, such that we can divide by two during the
 * interpolation without changing the result.
 *
 * @param c product ab with 2n-1 coefficients
 * @param a first multiplicand with n coefficients
 * @param b second multiplicand with n coefficients
 * @param n number of coefficients in a and b. needs to be divisble by 4
 * @param q modulus. needs to be co-prime to 3 and 5
 */
static void polymul_toom4(T* c, const T* a, const T* b, size_t n, T q){
    //TODO add rcs to the errors
    if(n%4 != 0) {
        printf("ERROR: toom4 currently only supports n divisible by 4. please pad your polynomial.");
        return;
    }
    if(q%3 == 0 || q%5 == 0) {
        printf("ERROR: toom4 currently only supports q co-prime to 3 and 5.");
        return;
    }

    if(q>(1<<(sizeof(T)*8-3))){
        printf("ERROR: toom3 currently requires q to be at least three bit smaller than the word size.");
        return;
    }


    // During the interpolation, we need to divide by 8, 4, 2, and 3, 5
    // In case q is co-prime to 2 and 3, we can simply multiply by the inverse

    // As q is often a power of two, the inverse of two does often not exist.
    // In that case there is a workaround: We do the small multiplications
    // modulu 8q instead of q, such that we can divide by two during the
    // interpolation without changing the result.

    T inv8 = zq_inverse(8, q);
    T inv4 = zq_inverse(4, q);
    T inv2 = zq_inverse(2, q);
    T qq = q;

    if(inv8 == 0){
        qq = 8*q;
    }

    T inv3 = zq_inverse(3, qq);
    T inv5 = zq_inverse(5, qq);
    size_t nsplit = n/4;

    // split a and b in 4 parts each, i.e.,
    // a = a0 + y a1 + y^2 a2 + y^3 a3
    // b = b0 + y b1 + y^2 b2 + y^3 b3
    const T *a0 = &a[0*nsplit];
    const T *a1 = &a[1*nsplit];
    const T *a2 = &a[2*nsplit];
    const T *a3 = &a[3*nsplit];

    const T *b0 = &b[0*nsplit];
    const T *b1 = &b[1*nsplit];
    const T *b2 = &b[2*nsplit];
    const T *b3 = &b[3*nsplit];

    // Evaluate at y = {0, \infty , 1, -1 , 2, -2, -3}
    const T *a_0 = a0; const T *a_inf = a3;
    T a_1[nsplit], a_m1[nsplit], a_2[nsplit], a_m2[nsplit], a_3[nsplit];
    const T *b_0 = b0; const T *b_inf = b3;
    T b_1[nsplit], b_m1[nsplit], b_2[nsplit], b_m2[nsplit], b_3[nsplit];

    T tmp0, tmp1;
    for(size_t i=0; i<nsplit; i++){
        // a_0  = a0
        // a_inf= a3
        // a_1  = a0 +   a1 +   a2 +    a3
        // a_m1 = a0 -   a1 +   a2 -    a3
        // a_2  = a0 + 2*a1 + 4*a2 +  8*a3
        // a_m2 = a0 - 2*a1 + 4*a2 -  8*a3
        // a_3  = a0 + 3*a1 + 9*a2 + 27*a3
        tmp0 = zq_mod(a0[i] + a2[i], qq);
        tmp1 = zq_mod(a1[i] + a3[i], qq);

        a_1[i]  = zq_mod(tmp0 + tmp1, qq);
        a_m1[i] = zq_mod(tmp0 + qq - tmp1, qq);

        tmp0 = zq_mod(a0[i] + ((T2)4*a2[i]), qq);
        tmp1 = zq_mod(zq_mod((T2)2*a1[i], qq) + zq_mod((T2)8*a3[i], qq), qq);

        a_2[i] = zq_mod(tmp0 + tmp1, qq);
        a_m2[i] = zq_mod(tmp0 + qq - tmp1, qq);

        a_3[i] = zq_mod(a0[i] + zq_mod((T2)3*a1[i], qq) + zq_mod((T2)9*a2[i], qq) + zq_mod((T2)27*a3[i], qq), qq);

        // b_0  = b0
        // b_inf= b3
        // b_1  = b0 +   b1 +   b2 +    b3
        // b_m1 = b0 -   b1 +   b2 -    b3
        // b_2  = b0 + 2*b1 + 4*b2 +  8*b3
        // b_m2 = b0 - 2*b1 + 4*b2 -  8*b3
        // b_3  = b0 + 3*b1 + 9*b2 + 27*b3

        tmp0 = zq_mod(b0[i] + b2[i], qq);
        tmp1 = zq_mod(b1[i] + b3[i], qq);

        b_1[i]  = zq_mod(tmp0 + tmp1, qq);
        b_m1[i] = zq_mod(tmp0 + qq - tmp1, qq);

        tmp0 = zq_mod(b0[i] + ((T2)4*b2[i]), qq);
        tmp1 = zq_mod(zq_mod((T2)2*b1[i], qq) + zq_mod((T2)8*b3[i], qq), qq);

        b_2[i] = zq_mod(tmp0 + tmp1, qq);
        b_m2[i] = zq_mod(tmp0 + qq - tmp1, qq);

        b_3[i] = zq_mod(b0[i] + zq_mod((T2)3*b1[i], qq) + zq_mod((T2)9*b2[i], qq) + zq_mod((T2)27*b3[i], qq), qq);
    }

    // do smaller multiplications (in this case using schoolbook; could also be
    // other Toom or Karatsuba)
    T c_0[2*nsplit-1], c_1[2*nsplit-1], c_m1[2*nsplit-1],
      c_2[2*nsplit-1], c_m2[2*nsplit-1], c_3[2*nsplit-1],
      c_inf[2*nsplit-1];

    poly_polymul_ref(c_0,  a_0, nsplit, b_0, nsplit, qq);
    poly_polymul_ref(c_1,  a_1, nsplit, b_1, nsplit, qq);
    poly_polymul_ref(c_m1, a_m1, nsplit, b_m1, nsplit, qq);
    poly_polymul_ref(c_2,  a_2, nsplit, b_2, nsplit, qq);
    poly_polymul_ref(c_m2, a_m2, nsplit, b_m2, nsplit, qq);
    poly_polymul_ref(c_3,  a_3, nsplit, b_3, nsplit, qq);
    poly_polymul_ref(c_inf, a_inf, nsplit, b_inf, nsplit, qq);

    // Interpolation of c = c_0 + v1 y^1 + v2 y^2 + v3 y^3 + v4 y^4 + v5 y^5 +
    //                      c_inf y^6

    for(size_t i=0; i<2*n-1; i++){
        c[i] = 0;
    }
    // need to use 16*q if no inverses for 2,4,8
    T2 qqq = (T2)q*16;
    T2 t0, t1, t2, v1, v2, v3, v4, v5;
    for(size_t i=0;i<2*nsplit-1;i++){
        if(inv8 == 0){
            // t0 = v2 +  v4
            t0 = zq_mod2((c_1[i] + c_m1[i])>>1, qqq);
            t0 = zq_mod2(t0 + qqq - c_0[i], qqq);
            t0 = zq_mod2(t0 + qqq - c_inf[i], qqq);
            // t1 = v2 + 4v4
            t1 = zq_mod2(c_2[i] + c_m2[i], qqq);
            t1 = zq_mod2(t1 + qqq - zq_mod2((T2)2*c_0[i], qqq), qqq);
            t1 = zq_mod2(t1 + qqq- zq_mod2((T2)128*c_inf[i], qqq), qqq);
            t1 = zq_mod2(t1>>3, qqq);

            v4 = zq_mod2((t1 + qqq - t0)*inv3, qqq);
            v2 = zq_mod2((t0 + qqq - v4), qqq);

            // t0 = v1 + v3 + v5
            t0 = zq_mod2(c_1[i] + qqq - c_m1[i], qqq);
            t0 = zq_mod2(t0 >> 1, qqq);
            // t1 = v3 + 5v5
            t1 = zq_mod2(c_2[i] + qqq - c_m2[i], qqq);
            t1 = zq_mod2(t1 >> 2, qqq);
            t1 = zq_mod2(t1 + qqq - t0, qqq);
            t1 = zq_mod2(t1 * inv3, qqq);
            // t2 = v1 + 9v3 + 81v5
            t2 = zq_mod2(c_3[i] + qqq - c_0[i], qqq);
            t2 = zq_mod2(t2 + qqq - zq_mod2((T2)9*v2, qqq), qqq);
            t2 = zq_mod2(t2 + qqq - zq_mod2((T2)81*v4, qqq), qqq);
            t2 = zq_mod2(t2 + qqq - zq_mod2((T2)729*c_inf[i], qqq), qqq);
            t2 = zq_mod2(t2*inv3, qqq);
            // t2 = 5v5
            t2 = zq_mod2(t2 + qqq - t0, qqq);
            t2 = zq_mod2(t2>>3, qqq);
            t2 = zq_mod2(t2 + qqq - t1, qqq);

            v5 = zq_mod2(t2*inv5, qqq);
            v3 = zq_mod2(t1 + qqq - t2, qqq);
            v1 = zq_mod2(t0 + qqq - v3, qqq);
            v1 = zq_mod2(v1 + qqq - v5, qqq);
        } else {
            // t0 = v2 +  v4
            t0 = zq_mod((c_1[i] + c_m1[i]) * inv2, qq);
            t0 = zq_mod(t0 + qq - c_0[i], qq);
            t0 = zq_mod(t0 + qq - c_inf[i], qq);
            // t1 = v2 + 4v4
            t1 = zq_mod(c_2[i] + c_m2[i], qq);
            t1 = zq_mod(t1 + 2*qq - 2*c_0[i], qq);
            t1 = zq_mod(t1 + qq - zq_mod((T2)128*c_inf[i], qq), qq);
            t1 = zq_mod(t1 * inv8, qq);

            v4 = zq_mod((t1 + qq - t0)*inv3, qq);
            v2 = zq_mod((t0 + qq - v4), qq);

            // t0 = v1 + v3 + v5
            t0 = zq_mod(c_1[i] + qq - c_m1[i], qq);
            t0 = zq_mod(t0 * inv2, qq);

            // t1 = v3 + 5v5
            t1 = zq_mod(c_2[i] + qq - c_m2[i], qq);
            t1 = zq_mod(t1 * inv4, qq);
            t1 = zq_mod(t1 + qq - t0, qq);
            t1 = zq_mod(t1 * inv3, qq);
            // t2 = v1 + 9v3 + 81v5
            t2 = zq_mod(c_3[i] + qq - c_0[i], qq);
            t2 = zq_mod(t2 + qq - zq_mod((T2)9*v2, qq), qq);
            t2 = zq_mod(t2 + qq - zq_mod((T2)81*v4, qq), qq);
            t2 = zq_mod(t2 + qq - zq_mod((T2)729*c_inf[i], qq), qq);
            t2 = zq_mod(t2*inv3, qq);
            // t2 = 5v5
            t2 = zq_mod(t2 + qq - t0, qq);
            t2 = zq_mod(t2*inv8, qq);
            t2 = zq_mod(t2 + qq - t1, qq);

            v5 = zq_mod(t2*inv5, qq);
            v3 = zq_mod(t1 + qq - t2, qq);
            v1 = zq_mod(t0 + qq - v3, qq);
            v1 = zq_mod(v1 + qq - v5, qq);
        }

        c[0*nsplit + i] = zq_mod(c[0*nsplit + i]+c_0[i], q);
        c[1*nsplit + i] = zq_mod(c[1*nsplit + i]+v1, q);
        c[2*nsplit + i] = zq_mod(c[2*nsplit + i]+v2, q);
        c[3*nsplit + i] = zq_mod(c[3*nsplit + i]+v3, q);
        c[4*nsplit + i] = zq_mod(c[4*nsplit + i]+v4, q);
        c[5*nsplit + i] = zq_mod(c[5*nsplit + i]+v5, q);
        c[6*nsplit + i] = zq_mod(c[6*nsplit + i]+c_inf[i], q);
    }
}

/**
 * @brief Random test of three-way Toom--Cook multiplication (`polymul_toom3`)
 *
 * @param n number of coefficients of input polynomials. needs to be divisble by 3
 * @param q modulus. needs to be co-prime to 3
 * @param printPoly flag for printing inputs and outputs
 * @return int  0 if test is successful, 1 otherwise
 */
static int testcase_toom3(size_t n, T q, int printPoly){
    int rc = 0;
    T a[n], b[n];
    T c_ref[2*n-1], c[2*n-1];
    printf("Testing Toom3 mul with q=%d, n=%ld\n", q, n);

    poly_random(a, n, q);
    poly_random(b, n, q);
    if(printPoly) poly_print("a", a, n);
    if(printPoly) poly_print("b", b, n);

    // compute reference product for comparison
    poly_polymul_ref(c_ref, a, n, b, n, q);
    if(printPoly) poly_print("c_ref", c_ref, 2*n-1);

    polymul_toom3(c, a, b, n, q);
    if(printPoly) poly_print("c", c, 2*n-1);

    if(poly_compare(c_ref, c, 2*n-1, q)) rc = 1;
    return rc;
}

/**
 * @brief Random test of three-way Toom--Cook multiplication (`polymul_toom4`)
 *
 * @param n number of coefficients of input polynomials. needs to be divisble by 4
 * @param q modulus needs to be co-prime to 3 and 5
 * @param printPoly flag for printing inputs and outputs
 * @return int  0 if test is successful, 1 otherwise
 */
static int testcase_toom4(size_t n, T q, int printPoly){
    int rc = 0;
    T a[n], b[n];
    T c_ref[2*n-1], c[2*n-1];
    printf("Testing Toom4 mul with q=%d, n=%ld\n", q, n);

    poly_random(a, n, q);
    poly_random(b, n, q);
    if(printPoly) poly_print("a", a, n);
    if(printPoly) poly_print("b", b, n);

    // compute reference product for comparison
    poly_polymul_ref(c_ref, a, n, b, n, q);
    if(printPoly) poly_print("c_ref", c_ref, 2*n-1);

    polymul_toom4(c, a, b, n, q);
    if(printPoly) poly_print("c", c, 2*n-1);

    if(poly_compare(c_ref, c, 2*n-1, q)) rc = 1;
    return rc;
}

int main(void){
    int rc = 0;

    // toom 3
    rc |= testcase_toom3(6, 17, 1);     // inverses of 2 and 3 exist
    rc |= testcase_toom3(6, 8, 1);      // inverse of 2 doesn't exist
    rc |= testcase_toom3(258, 1<<13, 0);// inverse of 2 doesn't exist
    rc |= testcase_toom3(258, 1<<15, 0);// inverse of 2 doesn't exist

    // toom 4
    rc |= testcase_toom4(8, 17, 1);     // inverses of 2,3,4,5,8 exist
    rc |= testcase_toom4(8, 1<<11, 1);  // inverses of 2,4,8 don't exist
    rc |= testcase_toom4(8, 1<<13, 1);  // inverses of 2,4,8 don't exist
    rc |= testcase_toom4(256, 17, 0);   // inverses of 2,3,4,5,8 exist
    rc |= testcase_toom4(256, 1<<11, 0);  // inverses of 2,4,8 don't exist
    rc |= testcase_toom4(256, 1<<12, 0);  // inverses of 2,4,8 don't exist
    // Saber parameters
    rc |= testcase_toom4(256, 1<<13, 0);// inverses of 2,4,8 don't exist
    rc |= testcase_toom4(256, 8191, 0); // inverses of 2,4,8 exist
    // TODO: get large moduli that are not a power-of-two to work
    //rc |= testcase_toom4(256, (1<<13)-16, 0);// inverses of 2,4,8 don't exist

    if(rc){
        printf("ERROR.\n");
    } else {
        printf("ALL GOOD.\n");
    }

    return rc;
}
