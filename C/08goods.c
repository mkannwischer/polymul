#include "common.h"
#include "poly.h"
#include "zq.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static void goodsPermutation(T *g, T *p, size_t p0, size_t p1){
    size_t i;
    for(i=0;i<p0*p1;i++){
        // g[i%p1].coeffs[i%p0] = p.coeffs[i]
        g[(i%p1)*p0 + i%p0] = p[i];
    }
}

static void goodsPermutationInverse(T *p, T *g, size_t p0, size_t p1){
    size_t i;
    for(i=0;i<p0*p1;i++){
        // p.coeffs[i] = g[i%p1].coeffs[i%p0]
        p[i] = g[(i%p1)*p0 + i%p0];
    }
}

static void goods(T *c, T *a, T *b, size_t p0, size_t p1, T q){
    size_t i,j;
    // find an p0-th root of unity for cyclic p0-NTT
    T root = zq_primitiveRootOfUnity(p0, q);

    // ONLINE COMPUTATION
    // compute Good's permutation
    T agood[p1][p0];
    T bgood[p1][p0];
    // Note: This trick can only be really fast if the Good's permutation
    // is done implicitly, i.e., as a part of the NTT.
    goodsPermutation(agood[0], a, p0, p1);
    goodsPermutation(bgood[0], b, p0, p1);

    // transform each element of a and b to NTT domain
    for(i=0;i<p1;i++){
        polymul_cyclic_ntt_forward_reference(agood[i], agood[i], q, p0, root);
        polymul_cyclic_ntt_forward_reference(bgood[i], bgood[i], q, p0, root);
    }

    // perform polynomial multiplication modulo x^p1-1
    // Alternatively, one could perform a p1-NTT here and then do a pointwise multiplication
    T cgood[p1][p0];
    for(i=0;i<p0;i++){
        T apoly[p1];
        T bpoly[p1];
        for(j=0; j<p1; j++){
            apoly[j] = agood[j][i];
            bpoly[j] = bgood[j][i];
        }

        T cpoly[2*p1-1];
        poly_polymul_ref(cpoly, apoly, p1, bpoly, p1, q);

        // reduce mod x^p1-1
        for(j=p1; j<2*p1-1; j++){
            cpoly[j-p1] = (cpoly[j-p1] + cpoly[j]) % q;
        }

        // write back to corresponding coefficent of c
        for(j=0;j<p1;j++){
            cgood[j][i] = cpoly[j];
        }
    }

    // perform inverse NTT
    for(i=0;i<p1;i++){
        polymul_cyclic_ntt_inverse_reference(cgood[i], cgood[i], q, p0, root);
    }

    goodsPermutationInverse(c, cgood[0], p0, p1);
}

static int testcase_goods(size_t p0, size_t p1, T q, int printPoly){
    int rc = 0;
    size_t n = p0*p1;
    T a[n], b[n];
    T c_ref[2*n-1], c[n];

    printf("Testing polynomial multiplication using Good's trick with n=%zu, p0=%zu, p1=%zu, q=%u\n", n, p0, p1, q);

    poly_random(a, n, q);
    poly_random(b, n, q);
    if(printPoly) poly_print("a=", a, n);
    if(printPoly) poly_print("b=", b, n);

    // compute reference product
    poly_polymul_ref(c_ref, a, n, b, n, q);
    // reduce mod x^n-1
    for(size_t i=n;i<2*n-1;i++){
        c_ref[i-n] = (c_ref[i-n] + c_ref[i]) % q;
    }
    if(printPoly) poly_print("a*b (ref)=", c_ref, n);
    goods(c, a, b, p0, p1, q);
    if(printPoly) poly_print("c=", c, n);
    if(poly_compare(c_ref, c, n, q)) rc = 1;
    return rc;
}

int main(void){
    int rc = 0;
    rc |= testcase_goods(8, 3, 17, 1);     // n=8*3=24, q=17
    rc |= testcase_goods(512, 3, 7681, 0); // n=512*3=1536, q=7681
    if(rc){
        printf("ERROR\n");
    } else {
        printf("ALL GOOD\n");
    }
    return rc;
}
