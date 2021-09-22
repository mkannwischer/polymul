/**
 * @file common.c
 * @brief Common code for pre-computation.
 *
 */
#include "common.h"

/**
 * @brief Bitreverse an array of length n inplace
 *
 * @param src array
 * @param n length of array
 */
void bitreverse(T *src, size_t n){
    for(size_t i = 0, j = 0; i < n; i++){
        if(i < j){
            src[i] += src[j];
            src[i] -= (src[j] = (src[i] - src[j]));
        }
        for(size_t k = n >> 1; (j ^= k) < k; k >>=1);
    }
}

/**
 * @brief Computes logarithm with arbitrary base
 *
 * @param x input
 * @param n base
 * @return unsigned int log_n(x)
 */
unsigned int log_base(unsigned int x, unsigned int n){
    double logd = log(x)/log(n);
    return (unsigned int) logd;
}
