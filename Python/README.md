
# The Poly datatype

We implement a `Poly` in [poly.py](./poly.py) class to represent all polynomials.
It consists of a list of `n` coefficients `coeffs` and optionally
a modulus `q` if the coefficient field/ring is modulo `q`.


A polynomial can be created as

```
a = Poly([1,2,3], q=17) # 1 + 2x + 3x^2
b = Poly([5,6,7], q=17) # 5 + 6x + 7x^2


c = Poly.zero(3, 17)    # all zero coefficients, n=3, q=17
d = Poly.random(3, 17)  # random coefficients in [0,q), n=3, q=17
```

We implement standard arithmetic
```
a = Poly([1,2,3], q=17)
b = Poly([5,6,7], q=17)

c = a + b
d = 17 * a
e = a - b
f = a * b  # results in a 2n-1 coefficient Poly
g = a >> 2


f.reduce() # reduction modulo q

c == a  # equality checks are checking coefficient vectors for equality
```


In addition, we implement a pointwise multiplication:
```
a = Poly([1,2,3], q=17)
b = Poly([5,6,7], q=17)

c = a.pointwise(b)  # coefficient-wise multiplication
```

# Common helper function

We implement a couple of common helper functions in [common.py](./common.py).
Most of them are required for the number-theoretic transforms.

Here is a list with brief descriptions. For more details see the documentation in [common.py](./common.py).
```
isPrime(n)                      # checks if a number is a prime
primitiveRootOfUnity(n, q)      # returns a primitive n-th root of unity modulo q

bitreverse(a)                   # bitreverses a list of length 2^k for some k
reverseBaseN(a, base)           # similiar as bitreverse, but for arbitrary base. a needs to be a list of length base^k for some k


# reference NTT implementations (n^2 complexity)
ntt_naive_cyclic(a, root)             # computes cyclic ntt of a; given n-th root of unity root
invntt_naive_cyclic(antt, root)       # computes inverse cyclic ntt of antt; given n-th root of unity root
ntt_naive_negacyclic(a, root)         # computes negacyclic ntt of a; given 2n-th root of unity root
invntt_naive_negacyclic(antt, root)   # computes inverse negacyclic ntt of anttl given 2n-th root of unity root
```