
# Framework

We represent polynomials by coefficient arrays.
The code is developed to support a generic datatype `T` for the coefficients,
which can be both `uint16_t` or `uint32_t`. Multiplications double the size to
`T2` which is either `uint32_t` or uint64_t`.
See [common.h](./common.h) for the type definitions.



## Poly helper functions ([poly.c](./poly.c))

For working with polynomials (represented by `T[]`), we implement some helper 
functions in [poly.c](./poly.c)):
```
void poly_print(char *polyname, T *a, size_t n);   // dumps a polynomial to stdout, e.g., "a= 16x^2 + 1x^1"
void poly_print2(T *a, size_t n);                  // dumps a polynomial to stdout, e.g., [0,1,16]

// performs a schoolbook multiplication computing the full product with 2*n-1 coeffs
void poly_polymul_ref(T *c, const T *a, size_t aN, const T *b, size_t bN, T q);

// samples a uniformly random polynomial with coefficients in [0, q)
void poly_random(T *a, size_t n, T q);

// compares two polynomials mod q
int poly_compare(T *a, T *b, size_t n, T q);

// slow reference implementations of cyclic NTTs and negacyclic NTTs
void polymul_cyclic_ntt_forward_reference(T *t, const T *a, T q, T n, T root);
void polymul_cyclic_ntt_inverse_reference(T *t, const T *a, T q, T n, T root);
void polymul_negacyclic_ntt_forward_reference(T *t, const T *a, T q, T n, T root);
void polymul_negacyclic_ntt_inverse_reference(T *t, const T *a, T q, T n, T root);
``

## Finite field helper functions ([fq.c](./fq.c))

We implement various helper functions for working in a finite field modulo q.
The common code can be found in [fq.c](./fq.c):

```
T zq_mod(T2 a, T q);                    // reduces a to [0, q) if q != 0; no-op if q = 0 (works only for small q).
T2 zq_mod2(T2 a, T2 q);                 // reduces a to [0, q) if q != 0; no-op if q = 0 (works for large q as well).
T zq_inverse(T a, T2 q);                // computes the multiplicative inverse of a mod q.
int zq_isPrime(T a);                    // returns 1 if a is a prime; 0 otherwise.
T zq_pow(T a, T2 e, T q);               // computes a^e mod q.
T zq_primitiveRootOfUnity(T n, T q);    // finds a primitive n-th root of unity modulo q.
```

Note that does are only used for precomputation. Hence, performance does not matter.

## Common helper functions ([common.c](./common.c))

There are a few helper functions which are not related to polynomials or finite
fields. Those are located in [common.c](./common.c):

```
void bitreverse(T *src, size_t n);                          // bitreverse an array of length n inplace
unsigned int log_base(unsigned int x, unsigned int n);      // compute logarithm of x in base n
```

Note that does are only used for precomputation. Hence, performance does not matter.