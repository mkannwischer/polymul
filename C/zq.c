#include "zq.h"


// naive implementation of modular reduction
// in the real life the implementation of this heavily depends on q
// in case q is a power of two, it is a nop
// in case it is a prime one would usually use Montgomery, Barrett, or
// specialized reductions for e.g., Solinas primes.
T zq_mod(T2 a, T q){
    // q=0: no reduction needed (e.g., because q=2^16)
    if(q==0) return a;
    return a%q;
}


T2 zq_mod2(T2 a, T2 q){
    if (q==0) return a;
    return a%q;
}

// naive implementation of modular inverse
// In a real implementation this would be pre-computed anyway
T zq_inverse(T a, T2 q){
    if(q==0) q=1<<16;
    for(T2 i=1; i<q; i++){
        T2 p = (T2) a*i;
        if(p%q == 1){
            return i;
        }
    }
    return 0;
}

// naive implementation for checking if a number is a prime
// In a real implementation this would be pre-computed anyway
int zq_isPrime(T a){
    for(T i=2; i<a;i++){
        if(a%i == 0) return 0;
    }
    return 1;
}


// naive variable time implementation of square-and-multiply
// In a real implementation this would be pre-computed anyway
T zq_pow(T a, T2 e, T q){
    if(e == 0) return 1;
    a = a % q;
    if(e == 1) return a;

    // pow(root, q-1, q) = 1 due to Fermat's little theorem
    e = e % (q-1);

    T t = zq_pow(a, e>>1, q);
    t = ((T2)t * t) % q;
    if (e % 2 == 0) {
        return t;
    } else {
        return (a * t) % q;
    }
}

// naive search of a primitive n-th root of unity mod q
// In a real implementation this would be pre-computed anyway
T zq_primitiveRootOfUnity(T n, T q){
    int isPrimitive;
    for(T i=2; i<q; i++){
        // check if this is an n-th root of unity.
        if(zq_pow(i, n, q) == 1){
            isPrimitive = 1;
            // i is an n-th root of unity, but it may not be primitive
            // check by making sure i^j != 1 mod q for 1<=j<i

            for(T j=2; j<n;j++){
                if(zq_pow(i, j, q) == 1){
                    isPrimitive = 0;
                }
            }

            if(isPrimitive){
                return i;
            }
        }
    }
    return 0;
}
