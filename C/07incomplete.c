#include "common.h"
#include "poly.h"
#include "zq.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static int precomp_ct_cyclic(T *twiddles, size_t n, T root, T q, size_t numLayers){
    size_t logn = log2(n);
    if(n != (1U<<logn)){
        printf("ERROR: cyclic fft requires n to be a power of two");
        return 1;
    }
    if(((q-1) % (1<<numLayers)) != 0){
        printf("ERROR: cyclic ntt requires q-1 divisible by 2**numLayers");
        return 1;
    }
    if(!zq_isPrime(q)){
        printf("ERROR: cyclic fft requires prime q.");
        return 1;
    }

    //powers = [pow(root, i, q) for i in range(2**numLayers//2)]
    T powers[(1<<(numLayers-1))];
    powers[0] = 1;
    for(size_t i=1;i<(1U<<(numLayers-1));i++){
        powers[i] = ((T2) powers[i-1]*root) % q;
    }
    bitreverse(powers, (1<<(numLayers-1)));

    for(size_t i = 0; i < numLayers; i++){
        for(size_t j = 0; j<(1U<<i); j++){
            *twiddles = powers[j];
            twiddles++;
        }
    }
    return 0;
}

static void precomp_basemul_cyclic(T *twiddles, T root, T q, size_t numLayers){
    //twiddles = [pow(root, i, q) for i in range(2**numLayers)]
    twiddles[0] = 1;
    for(size_t i=1;i<(1U<<(numLayers));i++){
        twiddles[i] = ((T2) twiddles[i-1]*root) % q;
    }
    bitreverse(twiddles, 1U<<(numLayers));
}

static int precomp_gs_cyclic(T *twiddles, size_t n, T root, T q, size_t numLayers){
    size_t logn = log2(n);
    if(n != (1U<<logn)){
        printf("ERROR: cyclic fft requires n to be a power of two");
        return 1;
    }
    if(((q-1) % (1<<numLayers)) != 0){
        printf("ERROR: cyclic ntt requires q-1 divisible by 2**numLayers");
        return 1;
    }
    if(!zq_isPrime(q)){
        printf("ERROR: cyclic fft requires prime q.");
        return 1;
    }

    //powers = [pow(root, -i, q) for i in range(2**numLayers//2)]
    T powers[(1<<(numLayers-1))];
    powers[0] = 1;
    T rootInverse = zq_inverse(root, q);
    for(size_t i=1;i< (1U<<(numLayers-1));i++){
        powers[i] = ((T2) powers[i-1]*rootInverse) % q;
    }
    bitreverse(powers,(1<<(numLayers-1)));

    for(size_t i = 0; i < numLayers; i++){
        for(size_t j = 0; j<(1U<<(numLayers-1-i)); j++){
            *twiddles = powers[j];
            twiddles++;
        }
    }
    return 0;
}

static int precomp_ct_negacyclic(T *twiddles, size_t n, T root, T q, size_t numLayers){
    size_t logn = log2(n);
    if(n != (1U<<logn)){
        printf("ERROR: cyclic fft requires n to be a power of two");
        return 1;
    }
    if(((q-1) % (1<<(numLayers + 1))) != 0){
        printf("ERROR: cyclic ntt requires q-1 divisible by 2**(numLayers+1)");
        return 1;
    }
    if(!zq_isPrime(q)){
        printf("ERROR: cyclic fft requires prime q.");
        return 1;
    }

    //powers = [pow(root, i, q) for i in range(2**numLayers//2)]
    T powers[(1<<numLayers)];
    powers[0] = 1;
    for(size_t i=1;i<(1U<<numLayers);i++){
        powers[i] = ((T2) powers[i-1]*root) % q;
    }
    bitreverse(powers, 1<<numLayers);

    for(size_t i = 0; i < (1U<<numLayers)-1; i++){
        twiddles[i] = powers[i+1];
    }
    return 0;
}

static void precomp_basemul_negacyclic(T *twiddles, T root, T q, size_t numLayers){
    //twiddles = [pow(root, 2*i+1, q) for i in range(2**numLayers)]
    twiddles[0] = root;
    T root2  = ((T2)root*root) % q;
    for(size_t i=1;i<(1U<<(numLayers));i++){
        twiddles[i] = ((T2) twiddles[i-1]*root2) % q;
    }
    bitreverse(twiddles, 1U<<(numLayers));
}

static int precomp_gs_negacyclic(T *twiddles, size_t n, T root, T q, size_t numLayers){
    size_t logn = log2(n);
    if(n != (1U<<logn)){
        printf("ERROR: cyclic fft requires n to be a power of two");
        return 1;
    }
    if(((q-1) % (1<<(numLayers+1))) != 0){
        printf("ERROR: cyclic ntt requires q-1 divisible by 2**(numLayers+1)");
        return 1;
    }
    if(!zq_isPrime(q)){
        printf("ERROR: cyclic fft requires prime q.");
        return 1;
    }

    //powers = [pow(root, -(i+1), q) for i in range(2**numLayers)]
    T powers[(1<<numLayers)];
    T rootInverse = zq_inverse(root, q);
    powers[0] = rootInverse;
    for(size_t i=1;i< 1U<<numLayers;i++){
        powers[i] = ((T2) powers[i-1]*rootInverse) % q;
    }
    bitreverse(powers, 1<<numLayers);
    for(size_t i=0;i<(1U<<numLayers)-1;i++){
        twiddles[i] = powers[i];
    }
    return 0;
}


static void ntt_ct(T *a, T *twiddles, size_t n, T q, size_t numLayers){
    size_t logn = log2(n);

    for(size_t i=0; i < numLayers; i++){
        size_t distance = 1U<< (logn - 1 -i);
        for(size_t j=0; j<(1U<<i); j++){
            T twiddle = *twiddles;
            twiddles++;
            // Note: in the cyclic case many of the twiddles are 1;
            // could optimize those multiplications away
            for(size_t k =0; k<distance; k++){
                size_t idx0 = 2*j*distance + k;
                size_t idx1 = idx0 + distance;
                T a0  = a[idx0];
                T a1  = ((T2) a[idx1] * twiddle) % q;
                a[idx0] = (a0 + a1) % q;
                a[idx1] = (a0 + q - a1) % q;
            }
        }
    }
}

static void invntt_gs(T *a, T *twiddles, size_t n, T q, size_t numLayers){
    size_t logn = log2(n);
    for(size_t i=logn-numLayers; i < logn; i++){
        size_t distance = 1<<i;
        for(size_t j=0; j<(1U<<(logn - 1 -i)); j++){
            T twiddle = *twiddles;
            twiddles++;
            // Note: in the cyclic case many of the twiddles are 1;
            // could optimize those multiplications away
            for(size_t k =0; k<distance; k++){
                size_t idx0 = 2*j*distance + k;
                size_t idx1 = idx0 + distance;
                T a0  = (a[idx0] + a[idx1]) % q;
                T a1  = (a[idx0] + q - a[idx1]) % q;
                a[idx0] = a0;
                a[idx1] = ((T2)a1*twiddle) % q;
            }
        }
    }

    // Note: Half of these multiplications can be merged into the last
    // layer of butterflies by pre-computing (twiddle*ninv)%q
    T ninv = zq_inverse(1U<<numLayers, q);
    for(size_t i=0;i<n;i++){
        a[i] = ((T2)a[i]*ninv)%q;
    }
}

static void polymul_ntt_ct_gs(T *c, T *a, T *b, T *twiddlesNtt, T *twiddlesInvNtt,
                              T *twiddlesBasemul, T n, T q, size_t numLayers){
    size_t logn = log2(n);
    size_t pointwiseDegree = 1<<(logn-numLayers);

    // NTT(a) and NTT(b)
    ntt_ct(a, twiddlesNtt, n, q, numLayers);
    ntt_ct(b, twiddlesNtt, n, q, numLayers);

    // basemul c = a \circ b
    T cp[2*pointwiseDegree-1];
    for(size_t i=0;i<n/pointwiseDegree;i++){
        poly_polymul_ref(cp, &a[i*pointwiseDegree], pointwiseDegree, &b[i*pointwiseDegree], pointwiseDegree, q);
        for(size_t j=0; j<pointwiseDegree; j++){
            c[i*pointwiseDegree+j] = cp[j];
        }
        for(size_t j=pointwiseDegree; j<2*pointwiseDegree-1; j++){
            c[i*pointwiseDegree+j-pointwiseDegree] = (c[i*pointwiseDegree+j-pointwiseDegree] + (T2)cp[j]*twiddlesBasemul[i]) % q;
        }
    }
    // inverse ntt NTT^-1(c)
    invntt_gs(c, twiddlesInvNtt, n, q, numLayers);
}

static int testcase_cyclic(size_t n, T q, size_t numLayers, int printPoly){
    int rc = 0;
    T a[n], b[n];
    T c_ref[2*n-1], c[n];

    T twiddlesNtt[n-1];
    T twiddlesInvNtt[n-1];
    T twiddlesBasemul[n>>1];

    // find a (1<<numLayers)-th root of unity
    T root = zq_primitiveRootOfUnity(1<<numLayers, q);
    printf("Testing polynomial multiplication using cyclic incomplete NTT (%zu layers) with n=%zu, q=%d\n", numLayers, n, q);

    // precompute twiddle factors
    rc |= precomp_ct_cyclic(twiddlesNtt, n, root, q, numLayers);
    rc |= precomp_gs_cyclic(twiddlesInvNtt, n, root, q, numLayers);
    precomp_basemul_cyclic(twiddlesBasemul, root, q, numLayers);

    poly_random(a, n, q);
    poly_random(b, n, q);
    if(printPoly) poly_print("a=", a, n);
    if(printPoly) poly_print("b=", b, n);

    // compute reference product
    poly_polymul_ref(c_ref, a, n, b, n, q);
    // reduce mod x^n -1
    for(size_t i=n;i<2*n-1;i++){
        c_ref[i-n] = (c_ref[i-n] + c_ref[i]) % q;
    }
    if(printPoly) poly_print("a*b (ref)=", c_ref, n);
    polymul_ntt_ct_gs(c, a, b, twiddlesNtt, twiddlesInvNtt, twiddlesBasemul, n, q, numLayers);
    if(printPoly) poly_print("c=", c, n);
    if(poly_compare(c_ref, c, n, q)) rc = 1;
    return rc;
}

static int testcase_negacyclic(size_t n, T q, size_t numLayers, int printPoly){
    int rc = 0;
    T a[n], b[n];
    T c_ref[2*n-1], c[n];

    T twiddlesNtt[n-1];
    T twiddlesInvNtt[n-1];
    T twiddlesBasemul[n>>1];
    // find a (1<<(numLayers+1))-th root of unity
    T root = zq_primitiveRootOfUnity(1<<(numLayers + 1), q);
    printf("Testing polynomial multiplication using negacyclic incomplete NTT (%zu layers) with n=%zu, q=%d\n", numLayers, n, q);

    // precompute twiddle factors
    rc |= precomp_ct_negacyclic(twiddlesNtt, n, root, q, numLayers);
    rc |= precomp_gs_negacyclic(twiddlesInvNtt, n, root, q, numLayers);
    precomp_basemul_negacyclic(twiddlesBasemul, root, q, numLayers);

    poly_random(a, n, q);
    poly_random(b, n, q);
    if(printPoly) poly_print("a=", a, n);
    if(printPoly) poly_print("b=", b, n);

    // compute reference product
    poly_polymul_ref(c_ref, a, n, b, n, q);
    // reduce mod x^n + 1
    for(size_t i=n;i<2*n-1;i++){
        c_ref[i-n] = (c_ref[i-n] + q -  c_ref[i]) % q;
    }
    if(printPoly) poly_print("a*b (ref)=", c_ref, n);
    polymul_ntt_ct_gs(c, a, b, twiddlesNtt, twiddlesInvNtt, twiddlesBasemul, n, q, numLayers);
    if(printPoly) poly_print("c=", c, n);
    if(poly_compare(c_ref, c, n, q)) rc = 1;
    return rc;
}

int main(void){
    int rc = 0;

    // test cyclic NTT (mod x^n-1)
    rc |= testcase_cyclic(8, 5, 2, 1);      // n=8, q=5, layers=2
    rc |= testcase_cyclic(8, 5, 1, 1);      // n=8, q=5, layers=1
    rc |= testcase_cyclic(256, 3329, 7, 0); // n=256, q=3329, layers=7
    rc |= testcase_cyclic(256, 3329, 6, 0); // n=256, q=3329, layers=6

    // test negacyclic NTT (mod x^n+1)
    rc |= testcase_negacyclic(8, 17, 2, 1);     // n=8, q=17, layers=2
    rc |= testcase_negacyclic(8, 17, 1, 1);     // n=8, q=17, layers=1
    rc |= testcase_negacyclic(256, 3329, 7, 0); // n=256, q=3329, layers=7 (Kyber)
    rc |= testcase_negacyclic(256, 3329, 6, 0); // n=256, q=3329, layers=6

    if(rc){
        printf("ERROR\n");
    } else {
        printf("ALL GOOD\n");
    }
    return rc;
}
