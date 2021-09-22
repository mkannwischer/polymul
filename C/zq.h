/**
 * @file zq.h
 * @brief Common code for working in the finite field Z/qZ (Z_q).
 *
 * Most code is only for pre-computation and not very efficient.
 */
#ifndef ZQ_H
#define ZQ_H

#include "common.h"


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
T zq_mod(T2 a, T q);

/**
 * @brief reduces a to [0, q) if q != 0; no-op if q=0
 *
 * works also for large q
 *
 * @param a element to be reduced
 * @param q modulus
 * @return T2 a mod q in [0, q)
 */
T2 zq_mod2(T2 a, T2 q);

/**
 * @brief Computes inverse a^-1 of a, such that a*a^-1 = 1 mod q
 *
 * naive implementation; only used in pre-computation anyway
 *
 * @param a element to be inverted
 * @param q modulus
 * @return T inverse
 */
T zq_inverse(T a, T2 q);

/**
 * @brief Checks if a is a prime
 *
 * naive implementation; only used in pre-computation anyway
 *
 * @param a element to be checked
 * @return int 1 is a is prime; 0 otherwise
 */
int zq_isPrime(T a);

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
T zq_pow(T a, T2 e, T q);

/**
 * @brief Finds an primitive n-th root of unity mod q
 *
 * Naive search; only used in pre-computation anyway.
 *
 * @param n n
 * @param q modulus
 * @return T root, s.t., root^k = 1 mod q with root^l != 1 mod q for all l < k.
 */
T zq_primitiveRootOfUnity(T n, T q);
#endif
