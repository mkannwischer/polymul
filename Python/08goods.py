#!/usr/bin/env python3
"""Section 2.2.8: Good's trick.

Illustrates Good's trick.
"""
from poly import Poly
import common
import sys
import math


def goodsPermutation(p, p0, p1, q):
    """Compute Good's permuation of polynomial p.

    Maps a polynomial in Zq[x]/(x^(p0p1)-1) to Zq[y]/(y^p1 − 1)[z]/(z^p0 − 1) by
    setting x = yz.
    Under the hood this converts a p0p1-coefficient polynomial into p1 polynomials
    of p0 coefficients. A coefficient with index i, will be permuted into
    the (i % p0)-th coeffcieint of the (i % p1)-th polynomial.

    Parameters
    ----------
    p : Poly
        polynomial to be permuted.
    p0 : int
        n=p0p1.
    p1 : int
        n=p1p1.
    q: int
        modulus.
    Returns
    ----------
    Poly
        permuted polynomial.
    """
    assert p.n == p0*p1
    # permute polynomial into p1 polynomials of p0 coefficients each
    g = [Poly.zero(p0, q) for _ in range(p1)]
    for i in range(p.n):
        g[i%p1].coeffs[i%p0] = p.coeffs[i]
    return g

def goodsPermutationInverse(g, p0, p1, q):
    """Compute inverse Good's permuation of polynomial g.

    Inverse of `goodsPermutation`.
    Maps Zq[y]/(y^p1 − 1)[z]/(z^p0 − 1) to  Zq[x]/(x^(p0p1)-1) with x=yz.
    Under the hood this converts p1 polynomials of p0 coeffcients into
    one polynomial of p0p1 coefficients.
    Given an index i0 (index of coefficient) and i1 (index of polynomial), the
    index in the new polynomial is i, such that,
    i0 = i mod p0 and i1 = i mod p1, i.e.,

    i = (p1^-1 mod p0) p1 i0 + (p0^-1 mod p1) p0 i1
    Parameters
    ----------
    g : Poly
        polynomial to be permuted.
    p0 : int
        n=p0p1.
    p1 : int
        n=p1p1.
    q: int
        modulus.
    Returns
    ----------
    Poly
        permuted polynomial.
    """
    assert len(g) == p1
    assert g[0].n == p0
    p = Poly.zero(p0*p1, q)
    for i in range(p.n):
        p.coeffs[i] = g[i%p1].coeffs[i%p0]
    return p

def goods(a, b, p0, p1, q):
    """Compute polynomial product ab using Good's trick.

    Full workflow illustrating Good's trick:
    1. Apply Good's permutation to a and b.
    2. Compute p0-NTT to each of the p1 polynomials of both a and b
    3. Perform a base multiplication of a and b modulo (y^p1 - 1); p0 base multiplications in total.
    4. Compute inverse p0-NTT to each of the p1 polynomial of the result.
    5. Apply inverse Good's permutation to result.

    Parameters
    ----------
    a : Poly
        First multiplicand with p0p1 coefficients.
    b : Poly
        Second multiplicand with p0p1 coefficients.
    p0 : int
        n=p0p1.
    p1 : int
        n=p1p1.
    q: int
        modulus.
    Returns
    ----------
    Poly
        Product ab.
    """
    assert a.n == p0*p1
    assert b.n == p0*p1
    assert math.gcd(p0, p1) == 1
    assert common.isPrime(p1)
    assert ((q-1) % (p0)) == 0

    # find a p0-th root of unity
    root = common.primitiveRootOfUnity(p0, q)

    # apply Good's permutation to both a and b
    a = goodsPermutation(a, p0, p1, q)
    b = goodsPermutation(b, p0, p1, q)

    # transform each element of a and b to NTT domain
    for i in range(p1):
       # of course we would use FFTs here, usually.
       a[i] = common.ntt_naive_cyclic(a[i], root)
       b[i] = common.ntt_naive_cyclic(b[i], root)

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
        c[i] = common.invntt_naive_cyclic(c[i], root)

    # perform inverse Good's permutation
    c = goodsPermutationInverse(c, p0, p1, q)
    return c


def testcase_goods(p0, p1, q, printPoly=True):
    """Random test of Good's trick.

    Parameters
    ----------
    p0 : int
        n=p0p1 number of coefficients in polynomials.
    p1 : int
        n=p0p1 number of coefficients in polynomials.
    q : int
        modulus.
    numLayers: list
        number of layers in the NTT. Needs to be <= log n.
    printPoly : boolean
        flag for printing inputs and outputs.
    Returns
    ----------
    int
        0 if test is successful, 1 otherwise.
    """
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