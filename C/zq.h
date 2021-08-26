#ifndef ZQ_H
#define ZQ_H

#include "common.h"



T zq_mod(T2 a, T q);
T2 zq_mod2(T2 a, T2 q);

// naive implementation of modular inverse
// In a real implementation this would be pre-computed anyway
T zq_inverse(T a, T2 q);

// naive implementation for checking if a number is a prime
// In a real implementation this would be pre-computed anyway
int zq_isPrime(T a);

// naive variable time implementation of square-and-multiply
// In a real implementation this would be pre-computed anyway
T zq_pow(T a, T2 e, T q);

// naive search of a primitive n-th root of unity mod q
// In a real implementation this would be pre-computed anyway
T zq_primitiveRootOfUnity(T n, T q);


#endif