#include "poly.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

// performs a Karatsuba multiplication computing the full product with 2*n-1 coeffs
//         Works for odd n.
//         Smaller polynomial multiplications are implemented using schoolbook.
//         For example, for n=4.
//         Let a0b0 = w; a1b1 = y; (a0+a1)(b0+b1)= z

//         We compute:
//           0    1    2    3    4    5    6
//          ---- ---- ----      ---- ---- ----
//         | w0 | w1 | w2 |    | y0 | y1 | y2 |
//          ---- ---- ---- ---- ---- ----  ----
//                 + | z0 | z1 | z2 |
//                    ---- ---- ----
//                 - | w0 | w1 | w2 |
//                    ---- ---- ----
//                 - | y0 | y1 | y2 |
//                    ---- ---- ----
static void polymul_karatsuba(T* c, const T* a, const T* b, size_t n, T q){
    size_t nhalf = n>>1;
    size_t i;
    // buffer to hold ((a_0 + a_1)(b_0 + b_1) - a_1 b_1 - a_0 b_0)
    T tmp[2*(n-nhalf)-1];

    // in case N is odd, a0,b0 will be the shorter limbs
    const T *a0 = &a[0];
    const T *a1 = &a[nhalf];
    const T *b0 = &b[0];
    const T *b1 = &b[nhalf];

    // temporary buffers to hold (a_0 + a_1) and (b_0 + b_1)
    // we can use the output buffer to reduce RAM footprint
    T *a0a1 = &c[0];       // N-nhalf coefficients
    T *b0b1 = &c[n-nhalf]; // N-nhalf coefficients

    for(i=0;i<n-nhalf;i++) {
        a0a1[i] = a1[i];
        b0b1[i] = b1[i];
    }
    for(i=0;i<nhalf;i++){
        a0a1[i] = (a0a1[i] + a0[i]) % q;
        b0b1[i] = (b0b1[i] + b0[i]) % q;
    }

    // compute ((a_0 + a_1)(b_0 + b_1) - a_1 b_1 - a_0 b_0)
    poly_polymul_ref(tmp, a0a1, n-nhalf, b0b1, n-nhalf, q);

    // compute a0*b0 and place it in c[0], ..., c[nhalf*2-2]
    poly_polymul_ref(&c[0], a0, nhalf, b0, nhalf, q);

    // compute a1*b1 and place it in c[2*nhalf], ..., c[2*nhalf+2*(N-nhalf)-2]
    poly_polymul_ref(&c[2*nhalf], a1, n-nhalf, b1, n-nhalf, q);

    // one coefficeint is not occupied by a0b0 + a1b1 Y
    c[2*nhalf-1] = 0;

    // subtract a0b0
    for(i=0;i<2*nhalf-1;i++){
        tmp[i] = (tmp[i] + q - c[i]) % q;
    }
    // subtract a1b1
    for(i=0;i<2*(n-nhalf)-1;i++){
        tmp[i] = (tmp[i] + q - c[2*nhalf+i]) % q;
    }

    // add it to the large buffer
    for(i=0;i<2*(n-nhalf)-1;i++){
        c[i+nhalf] = (c[i+nhalf] + tmp[i]) % q;
    }
}

// Similar as polymul_karatsuba, but the smaller polynomial multiplications
// are also using karatsuba down to a threadhold from which we switch to
// schoolbook.
// Assumes n is the same for both a and b.

static void polymul_karatsuba_recursive(T* c, const T* a, const T *b,
                                        size_t n, T q, size_t threshold){
    size_t nhalf = n>>1;
    size_t i;
    // buffer to hold ((a_0 + a_1)(b_0 + b_1) - a_1 b_1 - a_0 b_0)
    T tmp[2*(n-nhalf)-1];

    if(n <= threshold){
        poly_polymul_ref(c, a, n, b, n, q);
        return;
    }

    // if we want to support odd n, we need to implement karatsuba that allows
    // to have differently sized inputs
    // Note that this not implies, nhalf = n-nhalf
    if((n&1) == 1){
        printf("ERR: N must not be odd, but is %ld\n", n);
        return;
    }

    const T *a0 = &a[0];
    const T *a1 = &a[nhalf];
    const T *b0 = &b[0];
    const T *b1 = &b[nhalf];

    // temporary buffers to hold (a_0 + a_1) and (b_0 + b_1)
    // we can use the output buffer to reduce RAM footprint
    T *a0a1 = &c[0];       // nhalf coefficients
    T *b0b1 = &c[nhalf];   // nhalf coefficients

    for(i=0;i<nhalf;i++) {
        a0a1[i] = (a0[i]+a1[i]) % q;
        b0b1[i] = (b0[i]+b1[i]) % q;
    }

    // compute ((a_0 + a_1)(b_0 + b_1) - a_1 b_1 - a_0 b_0)
    polymul_karatsuba_recursive(tmp, a0a1, b0b1, nhalf, q, threshold);

    // compute a0*b0 and place it in c[0], ..., c[nhalf*2-2]
    polymul_karatsuba_recursive(&c[0], a0, b0, nhalf, q, threshold);

    // compute a1*b1 and place it in c[2*nhalf], ..., c[2*nhalf+2*(N-nhalf)-2]
    polymul_karatsuba_recursive(&c[2*nhalf], a1, b1, nhalf, q, threshold);

    // one coefficeint is not occupied by a0b0 + a1b1 Y
    c[2*nhalf-1] = 0;

    // subtract a0b0
    for(i=0;i<2*nhalf-1;i++){
        tmp[i] = (tmp[i] + q - c[i]) % q;
    }
    // subtract a1b1
    for(i=0;i<2*nhalf-1;i++){
        tmp[i] = (tmp[i] + q - c[2*nhalf+i]) % q;
    }

    // add it to the large buffer
    for(i=0;i<2*nhalf-1;i++){
        c[i+nhalf] = (c[i+nhalf] + tmp[i]) % q;
    }
}
/*
    The core observation of refined Karatsuba is that some additions are performed twice.
    Consider the example with n =4.
    Let a0b0 = w; a1b1 = y; (a0+a1)(b0+b1)= z
    Normal Karatsuba computes:
      0    1    2    3    4    5    6
     ---- ---- ----      ---- ---- ----
    | w0 | w1 | w2 |    | y0 | y1 | y2 |
     ---- ---- ---- ---- ---- ----  ----
            + | z0 | z1 | z2 |
               ---- ---- ----
            - | w0 | w1 | w2 |
               ---- ---- ----
            - | y0 | y1 | y2 |
               ---- ---- ----
    Here w2-y0 is computed twice (in column 2 and 4).
    For larger polynomials, this duplicate computation becomes significant (n//2-1)
    We can, thus, compute that part h once and save some additions.
    Consequently, refined Karatsuba looks like
      0    1    2    3    4    5    6
     ---- ---- ----      ---- ---- ----
    | w0 | w1 | h0 |    |-h0 | y1 | y2 |
     ---- ---- ---- ---- ---- ----  ----
            + | z0 | z1 | z2 |
               ---- ---- ----
            - | w0 | w1 |
               ---- ---- ----
            -      | y1 | y2 |
                    ---- ----
*/
static void polymul_refined_karatsuba(T* c, const T* a, const T* b, size_t n,
                                      T q){
    size_t nhalf = n>>1;
    size_t i;
    // buffer to hold ((a_0 + a_1)(b_0 + b_1) - a_1 b_1 - a_0 b_0)
    T tmp[2*(n-nhalf)-1];

    // in case N is odd, a0,b0 will be the shorter limbs
    const T *a0 = &a[0];
    const T *a1 = &a[nhalf];
    const T *b0 = &b[0];
    const T *b1 = &b[nhalf];

    // temporary buffers to hold (a_0 + a_1) and (b_0 + b_1)
    // we can use the output buffer to reduce RAM footprint
    T *a0a1 = &c[0];       // N-nhalf coefficients
    T *b0b1 = &c[n-nhalf]; // N-nhalf coefficients

    for(i=0;i<n-nhalf;i++) {
        a0a1[i] = a1[i];
        b0b1[i] = b1[i];
    }
    for(i=0;i<nhalf;i++){
        a0a1[i] = (a0a1[i] + a0[i]) % q;
        b0b1[i] = (b0b1[i] + b0[i]) % q;
    }

    // compute ((a_0 + a_1)(b_0 + b_1) - a_1 b_1 - a_0 b_0)
    poly_polymul_ref(tmp, a0a1, n-nhalf, b0b1, n-nhalf, q);

    // compute a0*b0 and place it in c[0], ..., c[nhalf*2-2]
    poly_polymul_ref(&c[0], a0, nhalf, b0, nhalf, q);

    // compute a1*b1 and place it in c[2*nhalf], ..., c[2*nhalf+2*(N-nhalf)-2]
    poly_polymul_ref(&c[2*nhalf], a1, n-nhalf, b1, n-nhalf, q);

    // one coefficeint is not occupied by a0b0 + a1b1 Y
    c[2*nhalf-1] = 0;


    // compute the common sum h and place it in the right place
    for(i=0;i<nhalf-1;i++){
        c[nhalf+i] = (c[nhalf+i] + q - c[2*nhalf+i]) % q;
        c[2*nhalf+i] = q-c[nhalf+i];
    }

    // subtract lower part a0b0
    for(i=0;i<2*nhalf-nhalf;i++){
        tmp[i] = (tmp[i] + q - c[i]) % q;
    }
    // subtract upper part of a1b1
    for(i=0;i<2*(n-nhalf)-nhalf;i++){
        tmp[i+nhalf-1] = (tmp[i+nhalf-1] + q - c[2*nhalf+nhalf-1+i]) % q;
    }

    // add it to the large buffer
    for(i=0;i<2*(n-nhalf)-1;i++){
        c[i+nhalf] = (c[i+nhalf] + tmp[i]) % q;
    }
}

static int testcase_karatsuba(size_t n, T q, int printPoly){
    int rc = 0;
    T a[n], b[n];
    T c_ref[2*n-1], c[2*n-1];
    printf("Testing normal Karatsuba with n=%ld, q=%d\n", n, q);

    poly_random(a, n, q);
    poly_random(b, n, q);
    if(printPoly) poly_print("a", a, n);
    if(printPoly) poly_print("b", b, n);
    // compute reference product for comparison
    poly_polymul_ref(c_ref, a, n, b, n, q);
    if(printPoly) poly_print("c_ref", c_ref, 2*n-1);

    polymul_karatsuba(c, a, b, n, q);
    if(printPoly) poly_print("c", c, 2*n-1);

    if(poly_compare(c_ref, c, 2*n-1, q)) rc = 1;
    return rc;
}

static int testcase_karatsuba_recursive(size_t n, T q, size_t t, int printPoly){
    int rc = 0;
    T a[n], b[n];
    T c_ref[2*n-1], c[2*n-1];
    printf("Testing recursive Karatsuba with n=%ld, q=%d\n", n, q);
    poly_random(a, n, q);
    poly_random(b, n, q);
    if(printPoly) poly_print("a", a, n);
    if(printPoly) poly_print("b", b, n);
    // compute reference product for comparison
    poly_polymul_ref(c_ref, a, n, b, n, q);
    if(printPoly) poly_print("c_ref", c_ref, 2*n-1);

    polymul_karatsuba_recursive(c, a, b, n, q, t);
    if(printPoly) poly_print("c", c, 2*n-1);

    if(poly_compare(c_ref, c, 2*n-1, q)) rc = 1;
    return rc;
}

static int testcase_refined_karatsuba(size_t n, T q, int printPoly){
    int rc = 0;
    T a[n], b[n];
    T c_ref[2*n-1], c[2*n-1];
    printf("Testing refined Karatsuba with n=%ld, q=%d\n", n, q);

    poly_random(a, n, q);
    poly_random(b, n, q);
    if(printPoly) poly_print("a", a, n);
    if(printPoly) poly_print("b", b, n);
    // compute reference product for comparison
    poly_polymul_ref(c_ref, a, n, b, n, q);
    if(printPoly) poly_print("c_ref", c_ref, 2*n-1);

    polymul_refined_karatsuba(c, a, b, n, q);
    if(printPoly) poly_print("c", c, 2*n-1);

    if(poly_compare(c_ref, c, 2*n-1, q)) rc = 1;
    return rc;
}

int main(void){
    int rc = 0;
    // plain Karatsuba tests
    rc |= testcase_karatsuba(5, 17, 1);
    rc |= testcase_karatsuba(256, 1<<13, 0);

    // recursive Karatsuba tests
    rc |= testcase_karatsuba_recursive(8, 17, 2, 1);
    rc |= testcase_karatsuba_recursive(256, 1<<13, 4, 0);

    // refined Karatsuba tests
    rc |= testcase_refined_karatsuba(4, 17, 1);
    rc |= testcase_refined_karatsuba(256, 1<<13, 0);

    if(rc){
        printf("ERROR.\n");
    } else {
        printf("ALL GOOD.\n");
    }

    return rc;
}
