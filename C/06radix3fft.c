#include "common.h"
#include "poly.h"
#include "zq.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static size_t reverseIdxBaseN(size_t idx, size_t ndigits, size_t base){
    if(ndigits == 1) {
        return idx;
    } else {
        size_t digit = (idx % base) * pow(base, ndigits-1);
        return digit + reverseIdxBaseN(idx/base, ndigits-1, base);
    }
}

static void reverseBaseN(T *twiddles, size_t num, size_t base){
    T tmp[num];
    memcpy(tmp, twiddles, sizeof(tmp));

    size_t ndigits = log_base(num, base);
    for(size_t j =0; j<num; j++){
        size_t reversedIdx = reverseIdxBaseN(j, ndigits, base);
        twiddles[reversedIdx]   = tmp[j];
    }
}

// this is a slightly weird reversing since it is reversing pairs of twiddles
static void reverseBaseN2(T *twiddles, size_t num, size_t base){
    T tmp[2*num];
    memcpy(tmp, twiddles, sizeof(tmp));

    size_t ndigits = log_base(num+2, base);
    for(size_t j =0; j<num; j++){
        size_t reversedIdx = reverseIdxBaseN(j, ndigits, base);
        twiddles[2*reversedIdx]   = tmp[2*j];
        twiddles[2*reversedIdx+1] = tmp[2*j+1];
    }
}

static int precomp_ntt_cyclic(T *twiddles, size_t n, T root, T q){
    size_t logn = log_base(n, 3);
    size_t numTwiddles = 2*pow(3, logn-1);
    if(n != zq_pow(3, logn, q)){
        printf("ERROR: cyclic fft requires n to be a power of 3");
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


    T powers[numTwiddles];
    for(size_t j =0; j<(numTwiddles)>>1; j++){
        powers[2*j]   = zq_pow(root, j, q);
        powers[2*j+1] = zq_pow(root, 2*j, q);
    }
    // reverse the powers
    reverseBaseN2(powers, (numTwiddles)>>1, 3);

    // arrange twiddles per layer
    for(size_t i = 0; i < logn; i++){
        for(size_t j = 0; j< pow(3, i); j++){
            *twiddles = powers[2*j];
            twiddles++;
            *twiddles = powers[2*j+1];
            twiddles++;
        }
    }

    return 0;
}

static int precomp_invntt_cyclic(T *twiddles, size_t n, T root, T q){
    size_t logn = log_base(n, 3);
    size_t numTwiddles = 2*pow(3, logn-1);

    if(n != zq_pow(3, logn, q)){
        printf("ERROR: cyclic fft requires n to be a power of 3");
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

    T powers[numTwiddles];

    T rootinv = zq_inverse(root, q);
    for(size_t j =0; j<(numTwiddles)>>1; j++){
        powers[2*j]   = zq_pow(rootinv, j, q);
        powers[2*j+1] = zq_pow(rootinv, 2*j, q);
    }
    // reverse the powers
    reverseBaseN2(powers, (numTwiddles)>>1, 3);
    // arrange twiddles per layer
    for(size_t i = 0; i < logn; i++){
        for(size_t j = 0; j< pow(3, logn-1-i); j++){
            *twiddles = powers[2*j];
            twiddles++;
            *twiddles = powers[2*j+1];
            twiddles++;
        }
    }
    return 0;
}


static void ntt(T *a, T* twiddles, size_t n, T root3, T q){
    T root3sq = ((T2)root3*root3) % q;
    size_t logn = log_base(n, 3);
    for(size_t i=0; i<logn;i++){
        size_t distance = pow(3, logn-1-i);
        for(size_t j=0; j< pow(3, i); j++){
            T c0 = twiddles[0]; // c
            T c1 = twiddles[1]; // c^2
            twiddles += 2;
            // a_i     = a_i +     c a_{i+n} +     c^2 a_{i+2n}
            // a_{i+n} = a_i +   w c a_{i+n} + w^2 c^2 a_{i+2n}
            // a_{i+2n}= a_i + w^2 c a_{i+n} + w   c^2 a_{i+2n}
            // twiddles are [c, c^2]
            for(size_t k=0; k<distance; k++){
                size_t idx0 = 3*j*distance + k;
                size_t idx1 = idx0+distance;
                size_t idx2 = idx1+distance;
                T a0 = a[idx0];
                T a1 = ((T2)a[idx1]*c0) % q;
                T a2 = ((T2)a[idx2]*c1) % q;

                a[idx0] = (a0 + a1 + a2) % q;
                T t0 = ((T2) a1*root3) % q;
                T t1 = ((T2) a2*root3sq) % q;
                a[idx1] = (a0 + t0 + t1) % q;
                t0 = ((T2) a1*root3sq) % q;
                t1 = ((T2) a2*root3) % q;
                a[idx2] = (a0 + t0 + t1) % q;
            }
        }
    }
}

static void invntt(T *a, T* twiddles, size_t n, T root3, T q){
    T root3sq = ((T2)root3*root3) % q;
    size_t logn = log_base(n, 3);
    for(size_t i=0; i<logn;i++){
        size_t distance = pow(3, i);
        for(size_t j=0; j< pow(3, logn-1-i); j++){
            T c0 = twiddles[0]; // 1/c
            T c1 = twiddles[1]; // 1/c^2
            twiddles += 2;
            // a_i     = 1/3     (a_i +     a_{i+n} +     a_{i+2n})
            // a_{i+n} = 1/(3c)  (a_i + w^2 a_{i+n} + w   a_{i+2n})
            // a_{i+2n}= 1/(3c^2)(a_i + w   a_{i_b} + w^2 a_{i+2n})

            // We delay the multiplyications by 1/3 til the end
            for(size_t k=0; k<distance; k++){
                size_t idx0 = 3*j*distance + k;
                size_t idx1 = idx0+distance;
                size_t idx2 = idx1+distance;
                T a0 = a[idx0];
                T a1 = a[idx1];
                T a2 = a[idx2];
                // a_i     = 1/3     (a_i +     a_{i+n} +     a_{i+2n})
                a[idx0] = (a0+a1+a2) % q;

                // a_{i+n} = 1/(3c)  (a_i + w^2 a_{i+n} + w   a_{i+2n})
                T t0 = ((T2)a1*root3sq) % q;
                T t1 = ((T2)a2*root3) % q;
                a[idx1] = (a0 + t0 + t1) % q;
                a[idx1] = ((T2)a[idx1]*c0) % q;

                // a_{i+2n}= 1/(3c^2)(a_i + w   a_{i_b} + w^2 a_{i+2n})
                t0 = ((T2)a1*root3) % q;
                t1 = ((T2)a2*root3sq) % q;
                a[idx2] = (a0 + t0 + t1) % q;
                a[idx2] = ((T2)a[idx2]*c1) % q;
            }
        }
    }

    // Note: 2/3 of these multiplications can be merged with the last layer
    // of butterflies.
    T ninv = zq_inverse(n, q);
    for(size_t i=0;i<n;i++){
        a[i] = ((T2)a[i]*ninv)%q;
    }
}
static void polymul_ntt(T *c, T *a, T *b, T *twiddlesNtt, T *twiddlesInvNtt, T root3, T n, T q){
    // NTT(a) and NTT(b)
    ntt(a, twiddlesNtt, n, root3, q);
    ntt(b, twiddlesNtt, n, root3, q);
    // pointwise product c = a \circ b
    for(size_t i=0;i<n;i++){
        c[i] = ((T2)a[i]*b[i])%q;
    }
    // inverse ntt NTT^-1(c)
    invntt(c, twiddlesInvNtt, n, root3, q);
}

static int testcase_cyclic(size_t n, T q, int printPoly){
    int rc = 0;
    T a[n], b[n], a2[n];
    T c_ref[2*n-1], c[n];
    T anttref[n];

    T twiddlesNtt[n-1];
    T twiddlesInvNtt[n-1];
    size_t logn = log_base(n, 3);

    // find an n-th root of unity
    T root = zq_primitiveRootOfUnity(n, q);
    // find the corresponding 3-rd root of unity
    T root3 = zq_pow(root, pow(3, logn-1), q);

    // precompute twiddles
    rc |= precomp_ntt_cyclic(twiddlesNtt, n, root, q);
    rc |= precomp_invntt_cyclic(twiddlesInvNtt, n, root, q);

    printf("Testing forward cyclic NTT with n=%zu, q=%d\n", n, q);
    poly_random(a, n, q);
    if(printPoly) poly_print("a", a, n);

    polymul_cyclic_ntt_forward_reference(anttref, a, q, n, root);

    // compute ntt inplafce
    ntt(a, twiddlesNtt, n, root3, q);
    // output of ntt is base-3 reversed; reverse anttref for comparison
    reverseBaseN(anttref, n, 3);

    if(printPoly) poly_print("antt=", a, n);
    if(printPoly) poly_print("anttref=", anttref, n);
    if(poly_compare(anttref, a, n, q)) rc = 1;

    printf("Testing inverse cyclic NTT with n=%zu, q=%d\n", n, q);
    poly_random(a, n, q);
    if(printPoly) poly_print("antt=", a, n);
    polymul_cyclic_ntt_inverse_reference(anttref, a, q, n, root);
    // GS needs inputs in bitreversed order
    reverseBaseN(a, n, 3);
    // compute invntt inplace
    invntt(a, twiddlesInvNtt, n, root3, q);
    if(printPoly) poly_print("a", a, n);
    if(printPoly) poly_print("aref", anttref, n);

    if(poly_compare(anttref, a, n, q)) rc = 1;

    printf("Testing inverse a == invntt(ntt(a)) cyclic with n=%zu, q=%d\n", n, q);
    poly_random(a, n, q);
    if(printPoly) poly_print("a", a, n);
    // copy to different buffer, so we can compare later
    memcpy(a2, a, sizeof(a));
    // forward ntt
    ntt(a2, twiddlesNtt, n, root3, q);
    if(printPoly) poly_print("antt=", a2, n);
    // inverse ntt
    invntt(a2, twiddlesInvNtt, n, root3, q);
    poly_print("a2=", a2, n);
    if(poly_compare(a2, a, n, q)) rc = 1;

    printf("Testing polynomial multiplication cyclic with n=%zu, q=%d\n", n, q);
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
    polymul_ntt(c, a, b, twiddlesNtt, twiddlesInvNtt, root3, n, q);
    if(printPoly) poly_print("c=", c, n);
    if(poly_compare(c_ref, c, n, q)) rc = 1;
    return rc;
}

int main(void){
    int rc = 0;
    rc |= testcase_cyclic(3, 7, 1);
    rc |= testcase_cyclic(9, 37, 1);
    rc |= testcase_cyclic(27, 109, 1);
    if(rc){
        printf("ERROR\n");
    } else {
        printf("ALL GOOD\n");
    }
    return rc;
}
