#include "common.h"

// TODO: document this function
void bitreverse(T *src, size_t n){
    for(size_t i = 0, j = 0; i < n; i++){
        if(i < j){
            src[i] += src[j];
            src[i] -= (src[j] = (src[i] - src[j]));
        }
        for(size_t k = n >> 1; (j ^= k) < k; k >>=1);
    }
}


unsigned int log_base(unsigned int x, unsigned int n){
    double logd = log(x)/log(n);
    return (unsigned int) logd;
}