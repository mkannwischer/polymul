#!/usr/bin/env python3
from poly import Poly
from common import primitiveRootOfUnity, ntt_naive_cyclic, invntt_naive_cyclic, isPrime
import sys
import math

# TODO: write documentation for each function

def goodsPermutation(p, p0, p1, q):
    assert p.n == p0*p1
    # permute polynomial into p1 polynomials of p0 coefficients each
    g = [Poly.zero(p0, q) for _ in range(p1)]
    for i in range(p.n):
        g[i%p1].coeffs[i%p0] = p.coeffs[i]
    return g

def goodsPermutationInverse(g, p0, p1, q):
    assert len(g) == p1
    assert g[0].n == p0
    p = Poly.zero(p0*p1, q)
    for i in range(p.n):
        p.coeffs[i] = g[i%p1].coeffs[i%p0]
    return p

def goods(a, b, p0, p1, q):
    assert a.n == p0*p1
    assert b.n == p0*p1
    assert math.gcd(p0, p1) == 1
    assert isPrime(p0)
    assert ((q-1) % (p0)) == 0

    # find a p0-th root of unity
    root = primitiveRootOfUnity(p0, q)

    # apply Good's permutation to both a and b
    a = goodsPermutation(a, p0, p1, q)
    b = goodsPermutation(b, p0, p1, q)

    # transform each element of a and b to NTT domain
    for i in range(p1):
       # of course we would use FFTs here, usually.
       a[i] = ntt_naive_cyclic(a[i], root)
       b[i] = ntt_naive_cyclic(b[i], root)

    # perform polynomial multiplication modulo x^p1-1
    # Alternatively, one could perform a p1-NTT here and then do a pointwise multiplication
    c = [Poly.zero(p0, q) for _ in range(p1)]
    for i in range(p0):
        a_poly = Poly([a[j].coeffs[i] for j in range(p1)], q=q)
        b_poly = Poly([b[j].coeffs[i] for j in range(p1)], q=q)

        c_poly = a_poly*b_poly
        # reduce mod x^p1-1
        for j in range(p1, c_poly.n):
            c_poly.coeffs[j-p1] += c_poly.coeffs[j]
        c_poly.reduce()

        # write back to corresponding coeffients in c
        for j in range(p1):
            c[j].coeffs[i] = c_poly.coeffs[j]

    # perform inverse NTT
    for i in range(p1):
        c[i] = invntt_naive_cyclic(c[i], root)

    # perform inverse Good's permutation
    c = goodsPermutationInverse(c, p0, p1, q)
    return c


def testcase_goods(p0, p1, q, printPoly=True):
    rc = 0
    n  = p0*p1

    print(f"Testing polynomial multiplication using Good's trick with n={n}, p0={p0}, p1={p1}, q={q}")
    a = Poly.random(n, q)
    b = Poly.random(n, q)
    if printPoly: print("a=", a)
    if printPoly: print("b=", b)

    # compute reference product using schoolbook for comparison
    c_ref = a*b
    # reduce mod x^n - 1
    c_red = Poly(c_ref.coeffs[:n], q)
    for i in range(n, c_ref.n):
        c_red.coeffs[i-n] += c_ref.coeffs[i]
    c_red.reduce()

    if printPoly: print("a*b (ref)=", c_red)
    c = goods(a, b, p0, p1, q)

    if printPoly: print("c (Good's NTT)=", c)
    print(f"equal: {c == c_red}")
    if c != c_red:
        rc = 1
    return rc

if __name__ == "__main__":
    rc = 0

    rc |= testcase_goods(p0=8, p1=3, q=17)
    rc |= testcase_goods(p0=512, p1=3, q=7681, printPoly=False)

    if rc != 0:
        print("TEST FAILED.")
        sys.exit(1)
    print("ALL GOOD.")