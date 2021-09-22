/**
 * @file zq.c
 * @brief Common code for working in the finite field Z/qZ (Z_q).
 *
 * Most code is only for pre-computation and not very efficient.
 */

#include "zq.h"

/**
 * @brief reduces a to [0, q) if q != 0; no-op if q=0
 *
 * works only for small q.
 * In the real life the implementation of this heavily depends on q.
 * In case q is a power of two, it is a nop.
 * In case it is a prime one would usually use Montgomery, Barrett, or
 * specialized reductions for e.g., Solinas primes.
 *
 * @param a element to be reduced
 * @param q modulus
 * @return T a mod q in [0, q)
 */
T zq_mod(T2 a, T q){
    // q=0: no reduction needed (e.g., because q=2^16)
    if(q==0) return a;
    return a%q;
}

/**
 * @brief reduces a to [0, q) if q != 0; no-op if q=0
 *
 * works also for large q
 *
 * @param a element to be reduced
 * @param q modulus
 * @return T2 a mod q in [0, q)
 */
T2 zq_mod2(T2 a, T2 q){
    if (q==0) return a;
    return a%q;
}

/**
 * @brief Computes inverse a^-1 of a, such that a*a^-1 = 1 mod q
 *
 * naive implementation; only used in pre-computation anyway
 *
 * @param a element to be inverted
 * @param q modulus
 * @return T inverse
 */
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

/**
 * @brief Checks if a is a prime
 *
 * naive implementation; only used in pre-computation anyway
 *
 * @param a element to be checked
 * @return int 1 is a is prime; 0 otherwise
 */
int zq_isPrime(T a){
    for(T i=2; i<a;i++){
        if(a%i == 0) return 0;
    }
    return 1;
}

/**
 * @brief Computes exponentation a^e mod q
 *
 * Naive variable time implementation of square-and-multiply
 * Only used in pre-computation anyway.
 *
 * @param a base
 * @param e exponent
 * @param q moduluis
 * @return T a^e mod q
 */
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

/**
 * @brief Finds an primitive n-th root of unity mod q
 *
 * Naive search; only used in pre-computation anyway.
 *
 * @param n n
 * @param q modulus
 * @return T root, s.t., root^k = 1 mod q with root^l != 1 mod q for all l < k.
 */
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
