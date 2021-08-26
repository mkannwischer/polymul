#include "common.h"
#include "poly.h"
#include "zq.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>


static int precomp_ct_cyclic(T *twiddles, size_t n, T root, T q){
    // Note: Some twiddles could be reused; we need n//2 twiddles for cyclic
    // NTTs, and in layer i we use the first 2^i of them
    // However, in this function we pre-compute them for each layer separately,
    // which allows us to use the same code for a cyclic and a negacyclic NTT.
    // In a real implementation, one would of course not do this.
    size_t logn = log2(n);
    if(n != (1U<<logn)){
        printf("ERROR: cyclic fft requires n to be a power of two");
        return 1;
    }
    if(((q-1) % (n)) != 0){
        printf("ERROR: cyclic ntt requires q-1 divisible by n");
        return 1;
    }
    if(!zq_isPrime(q)){
        printf("ERROR: cyclic fft requires prime q.");
        return 1;
    }

    //powers = [pow(root, i, q) for i in range(n//2)]
    T powers[n>>1];
    powers[0] = 1;
    for(size_t i=1;i< n>>1;i++){
        powers[i] = ((T2) powers[i-1]*root) % q;
    }
    bitreverse(powers, n>>1);

    for(size_t i = 0; i < logn; i++){
        for(size_t j = 0; j<(1U<<i); j++){
            *twiddles = powers[j];
            twiddles++;
        }
    }
    return 0;
}

static int precomp_ct_negacyclic(T *twiddles, size_t n, T root, T q){
    size_t logn = log2(n);
    if(n != (1U<<logn)){
        printf("ERROR: negacyclic fft requires n to be a power of two");
        return 1;
    }
    if(((q-1) % (2*n)) != 0){
        printf("ERROR: negacyclic fft requires q-1 divisible by 2*n");
        return 1;
    }
    if(!zq_isPrime(q)){
        printf("ERROR: negacyclic fft requires prime q.");
        return 1;
    }

    T powers[n];
    powers[0] = 1;
    for(size_t i=1;i<n;i++){
        powers[i] = ((T2) powers[i-1]*root) % q;
    }
    bitreverse(powers, n);

    for(size_t i=0; i<n-1;i++){
        twiddles[i] = powers[i+1];
    }
    return 0;
}

static void ntt_ct(T *a, T *twiddles, size_t n, T q){
    size_t logn = log2(n);

    for(size_t i=0; i < logn; i++){
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


static int precomp_gs_cyclic(T *twiddles, size_t n, T root, T q){
    // Note: some twiddles could be reused; we need n//2 twiddles for cyclic
    // NTTs, and in layer i we use the first 2^(log_2(n)-1-i) of them
    // However, in this function we pre-compute them for each layer separately,
    // which allows us to use the same code for a cyclic and a negacyclic invNTT.
    // In a real implementation, one would of course not do this.
    size_t logn = log2(n);
    if(n != (1U<<logn)){
        printf("ERROR: cyclic fft requires n to be a power of two");
        return 1;
    }
    if(((q-1) % (n)) != 0){
        printf("ERROR: cyclic ntt requires q-1 divisible by n");
        return 1;
    }
    if(!zq_isPrime(q)){
        printf("ERROR: cyclic fft requires prime q.");
        return 1;
    }

    //powers = [pow(root, -i, q) for i in range(n//2)]
    T powers[n>>1];
    powers[0] = 1;
    T rootinv = zq_inverse(root, q);
    for(size_t i=1;i< n>>1;i++){
        powers[i] = ((T2) powers[i-1]*rootinv) % q;
    }
    bitreverse(powers, n>>1);

    for(size_t i = 0; i < logn; i++){
        for(size_t j = 0; j<(1U<<(logn-1-i)); j++){
            *twiddles = powers[j];
            twiddles++;
        }
    }
    return 0;
}


static int precomp_gs_negacyclic(T *twiddles, size_t n, T root, T q){
    size_t logn = log2(n);
    if(n != (1U<<logn)){
        printf("ERROR: negacyclic fft requires n to be a power of two");
        return 1;
    }
    if(((q-1) % (2*n)) != 0){
        printf("ERROR: negacyclic ntt requires q-1 divisible by 2n");
        return 1;
    }
    if(!zq_isPrime(q)){
        printf("ERROR: negacyclic fft requires prime q.");
        return 1;
    }
    //powers = [pow(root, -(i+1), q) for i in range(n)]
    T powers[n];
    T rootInverse = zq_inverse(root, q);
    powers[0] = rootInverse;
    for(size_t i=1;i<n;i++){
        powers[i] = ((T2)powers[i-1]*rootInverse) % q;
    }
    bitreverse(powers, n);
    for(size_t i=0;i<n-1;i++){
        twiddles[i] = powers[i];
    }
    return 0;
}



static void invntt_gs(T *a, T *twiddles, size_t n, T q){
    size_t logn = log2(n);
    for(size_t i=0; i < logn; i++){
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
    T ninv = zq_inverse(n, q);
    for(size_t i=0;i<n;i++){
        a[i] = ((T2)a[i]*ninv)%q;
    }
}

static void polymul_ntt_ct_gs(T *c, T *a, T *b, T *twiddlesNtt, T *twiddlesInvNtt, T n, T q){
    // NTT(a) and NTT(b)
    ntt_ct(a, twiddlesNtt, n, q);
    ntt_ct(b, twiddlesNtt, n, q);
    // pointwise product c = a \circ b
    for(size_t i=0;i<n;i++){
        c[i] = ((T2)a[i]*b[i])%q;
    }
    // inverse ntt NTT^-1(c)
    invntt_gs(c, twiddlesInvNtt, n, q);
}

static int testcase_cyclic(size_t n, T q, int printPoly){
    int rc = 0;
    T a[n], a2[n], b[n];
    T c_ref[2*n-1], c[n];
    T anttref[n];

    T twiddlesNtt[n-1];
    T twiddlesInvNtt[n-1];

    // find a n-th root of unity
    T root = zq_primitiveRootOfUnity(n, q);

    // precompute twiddle factors
    rc |= precomp_ct_cyclic(twiddlesNtt, n, root, q);
    rc |= precomp_gs_cyclic(twiddlesInvNtt, n, root, q);

    printf("Testing forward cyclic NTT using CT butterflies with n=%zu, q=%d\n", n, q);

    poly_random(a, n, q);
    if(printPoly) poly_print("a", a, n);

    polymul_cyclic_ntt_forward_reference(anttref, a, q, n, root);
    // compute ntt inplace
    ntt_ct(a, twiddlesNtt, n, q);

    // output of fft is bitreversed; bitreverse anttref for comparison
    bitreverse(anttref, n);

    if(printPoly) poly_print("antt\t", a, n);
    if(printPoly) poly_print("anttref", anttref, n);

    if(poly_compare(anttref, a, n, q)) rc = 1;

    printf("Testing inverse cyclic NTT using GS butterflies with n=%zu, q=%d\n", n, q);
    poly_random(a, n, q);
    if(printPoly) poly_print("antt", a, n);

    polymul_cyclic_ntt_inverse_reference(anttref, a, q, n, root);
    //GS needs inputs in bitreversed order
    bitreverse(a, n);
    // compute invntt inplace
    invntt_gs(a, twiddlesInvNtt, n, q);
    if(printPoly) poly_print("a", a, n);
    if(printPoly) poly_print("aref", anttref, n);
    if(poly_compare(anttref, a, n, q)) rc = 1;

    printf("Testing inverse a == invntt(ntt(a)) cyclic with n=%zu, q=%d\n", n, q);
    poly_random(a, n, q);
    if(printPoly) poly_print("a", a, n);
    // copy to different buffer, so we can compare later
    memcpy(a2, a, sizeof(a));
    // forward ntt
    ntt_ct(a2, twiddlesNtt, n, q);
    if(printPoly) poly_print("antt", a2, n);
    // inverse ntt
    invntt_gs(a2, twiddlesInvNtt, n, q);
    if(printPoly) poly_print("a2", a2, n);
    if(poly_compare(a2, a, n, q)) rc = 1;

    printf("Testing polynomial multiplication cyclic with n=%zu, q=%d\n", n, q);
    poly_random(a, n, q);
    poly_random(b, n, q);
    if(printPoly) poly_print("a", a, n);
    if(printPoly) poly_print("b", b, n);

    // compute reference product
    poly_polymul_ref(c_ref, a, n, b, n, q);
    // reduce mod x^n-1
    for(size_t i=n;i<2*n-1;i++){
        c_ref[i-n] = (c_ref[i-n] + c_ref[i]) % q;
    }
    if(printPoly) poly_print("c_ref", c_ref, n);
    polymul_ntt_ct_gs(c, a, b, twiddlesNtt, twiddlesInvNtt, n, q);

    if(printPoly) poly_print("c", c, n);
    if(poly_compare(c_ref, c, n, q)) rc = 1;

    return rc;
}

static int testcase_negacyclic(size_t n, T q, int printPoly){
    int rc = 0;
    T a[n], a2[n], b[n];
    T c_ref[2*n-1], c[n];
    T anttref[n];

    T twiddlesNtt[n-1];
    T twiddlesInvNtt[n-1];

    // find a 2n-th root of unity
    T root = zq_primitiveRootOfUnity(2*n, q);

    // precompute twiddle factors
    rc |= precomp_ct_negacyclic(twiddlesNtt, n, root, q);
    rc |= precomp_gs_negacyclic(twiddlesInvNtt, n, root, q);

    printf("Testing forward negacyclic NTT using CT butterflies with n=%zu, q=%d\n", n, q);

    poly_random(a, n, q);
    if(printPoly) poly_print("a", a, n);

    polymul_negacyclic_ntt_forward_reference(anttref, a, q, n, root);
    // compute ntt inplace
    ntt_ct(a, twiddlesNtt, n, q);

    // output of fft is bitreversed; bitreverse anttref for comparison
    bitreverse(anttref, n);

    if(printPoly) poly_print("antt\t", a, n);
    if(printPoly) poly_print("anttref", anttref, n);

    if(poly_compare(anttref, a, n, q)) rc = 1;

    printf("Testing inverse negacyclic NTT using GS butterflies with n=%zu, q=%d\n", n, q);
    poly_random(a, n, q);
    if(printPoly) poly_print("antt", a, n);

    polymul_negacyclic_ntt_inverse_reference(anttref, a, q, n, root);
    //GS needs inputs in bitreversed order
    bitreverse(a, n);
    // compute invntt inplace
    invntt_gs(a, twiddlesInvNtt, n, q);
    if(printPoly) poly_print("a", a, n);
    if(printPoly) poly_print("aref", anttref, n);
    if(poly_compare(anttref, a, n, q)) rc = 1;

    printf("Testing inverse a == invntt(ntt(a)) negacyclic with n=%zu, q=%d\n", n, q);
    poly_random(a, n, q);
    if(printPoly) poly_print("a", a, n);
    // copy to different buffer, so we can compare later
    memcpy(a2, a, sizeof(a));
    // forward ntt
    ntt_ct(a2, twiddlesNtt, n, q);
    if(printPoly) poly_print("antt", a2, n);
    // inverse ntt
    invntt_gs(a2, twiddlesInvNtt, n, q);
    if(printPoly) poly_print("a2", a2, n);
    if(poly_compare(a2, a, n, q)) rc = 1;

    printf("Testing polynomial multiplication negacyclic with n=%zu, q=%d\n", n, q);
    poly_random(a, n, q);
    poly_random(b, n, q);
    if(printPoly) poly_print("a", a, n);
    if(printPoly) poly_print("b", b, n);

    // compute reference product
    poly_polymul_ref(c_ref, a, n, b, n, q);
    // reduce mod x^n+1
    for(size_t i=n;i<2*n-1;i++){
        c_ref[i-n] = (c_ref[i-n] + q - c_ref[i]) % q;
    }
    if(printPoly) poly_print("c_ref", c_ref, n);
    polymul_ntt_ct_gs(c, a, b, twiddlesNtt, twiddlesInvNtt, n, q);

    if(printPoly) poly_print("c", c, n);
    if(poly_compare(c_ref, c, n, q)) rc = 1;

    return rc;
}

int main(void){
    int rc = 0;

    // test cyclic NTT (mod x^n-1)
    rc |= testcase_cyclic(8, 17, 1);
    rc |= testcase_cyclic(256, 3329, 0);

    // test negacyclic NTT (mod x^n+1)
    rc |= testcase_negacyclic(8, 17, 1);
    rc |= testcase_negacyclic(256, 7681, 0);

    if(rc){
        printf("ERROR\n");
    } else {
        printf("ALL GOOD\n");
    }
    return rc;
}
