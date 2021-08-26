#!/usr/bin/env python3
from poly import Poly
from common import isPrime, primitiveRootOfUnity
import sys

# TODO: write documentation for each function

def polymul_ntt(a, b):
    assert a.n == b.n
    assert a.q == b.q
    n = a.n
    q = a.q
    assert isPrime(q)
    assert ((q-1) % (2*n-1)) == 0

    # For implementing a general purpose NTT (not cyclic, negacyclic), we need
    # to evaluate a and b at 2n-1 roots of unity.

    def ntt_naive(p, root):
        pntt = Poly.zero(2*p.n-1, q)
        for i in range(2*p.n-1):
            for j in range(p.n):
                pntt.coeffs[i] = (pntt.coeffs[i] +  p.coeffs[j]*pow(root, i*j, q)) % q
        return pntt

    def invntt_naive(pntt, root):
        q = pntt.q
        p = Poly.zero(pntt.n, q)
        for i in range(pntt.n):
            for j in range(pntt.n):
                p.coeffs[i] = (p.coeffs[i] + pntt.coeffs[j]*pow(root, -i*j, q)) % q

        ninv = pow(p.n, -1, q)
        for i in range(p.n):
            p.coeffs[i] = (p.coeffs[i] * ninv) % q
        return p

    # we need a 2n-th root of unity
    root = primitiveRootOfUnity(2*n-1, q)

    # transform a and b to NTT domain
    antt = ntt_naive(a, root)
    bntt = ntt_naive(b, root)

    # point-wise multiplication
    cntt = antt.pointwise(bntt)

    # trnasform the result back to normal domain
    c = invntt_naive(cntt, root)

    return c


def polymul_cyclic_ntt(a, b):
    assert a.n == b.n
    assert a.q == b.q
    n = a.n
    q = a.q
    assert isPrime(q)
    assert ((q-1) % (n)) == 0


    # we need a n-th root of unity
    root = primitiveRootOfUnity(n, q)

    # For implementing a cyclic NTT we need to evaluate a and b at n roots of unity.

    def ntt_naive(p, root):
        q = p.q
        pntt = Poly.zero(p.n, q)
        for i in range(p.n):
            for j in range(p.n):
                pntt.coeffs[i] = (pntt.coeffs[i] +  p.coeffs[j]*pow(root, i*j, q)) % q
        return pntt

    def invntt_naive(pntt, root):
        q = pntt.q
        p = Poly.zero(pntt.n, q)
        for i in range(pntt.n):
            for j in range(pntt.n):
                p.coeffs[i] = (p.coeffs[i] + pntt.coeffs[j]*pow(root, -i*j, q)) % q
        ninv = pow(p.n, -1, q)
        for i in range(p.n):
            p.coeffs[i] = (p.coeffs[i] * ninv) % q
        return p

    # transform a and b to NTT domain
    antt = ntt_naive(a, root)
    bntt = ntt_naive(b, root)

    # point-wise multiplication
    cntt = antt.pointwise(bntt)

    # trnasform the result back to normal domain
    c = invntt_naive(cntt, root)

    return c

def polymul_negacyclic_ntt(a, b):
    assert a.n == b.n
    assert a.q == b.q
    n = a.n
    q = a.q
    assert isPrime(q)
    # We need a 2n-th root of unity rather than an n-th root of unity
    assert ((q-1) % (2*n)) == 0

    # we need a n-th root of unity
    root = primitiveRootOfUnity(2*n, q)

    def ntt_naive(p, root):
        q = p.q
        pntt = Poly.zero(p.n, q)
        for i in range(p.n):
            for j in range(p.n):
                pntt.coeffs[i] = (pntt.coeffs[i] +  p.coeffs[j]*pow(root, j, q)*pow(root, 2*i*j, q)) % q
        return pntt

    def invntt_naive(pntt, root):
        q = pntt.q
        p = Poly.zero(pntt.n, q)
        for i in range(pntt.n):
            for j in range(pntt.n):
                p.coeffs[i] = (p.coeffs[i] + pntt.coeffs[j]*pow(root, -2*i*j, q)) % q
        ninv = pow(p.n, -1, q)
        for i in range(p.n):
            p.coeffs[i] = (p.coeffs[i] * ninv * pow(root, -i, q)) % q
        return p

    # transform a and b to NTT domain
    antt = ntt_naive(a, root)
    bntt = ntt_naive(b, root)

    # point-wise multiplication
    cntt = antt.pointwise(bntt)

    # transform the result back to normal domain
    c = invntt_naive(cntt, root)

    return c


def testcase_ntt(n, q, printPoly=True):
    print(f"Testing naive NTT multiplication with n={n}, q={q}")
    a = Poly.random(n, q)
    b = Poly.random(n, q)
    if printPoly: print("a=", a)
    if printPoly: print("b=", b)

    # compute reference product using schoolbook for comparison
    c_ref = a*b
    if printPoly: print("a*b (ref)=", c_ref)
    c = polymul_ntt(a, b)
    if printPoly: print("a*b (NTT)=", c)
    print(f"equal: {c == c_ref}")
    if c == c_ref:
        return 0
    else:
        return 1

def testcase_cyclic_ntt(n, q, printPoly=True):
    print(f"Testing naive cyclic NTT multiplication (i.e., mod x^n - 1) with q={q}, n={n}")
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
    c = polymul_cyclic_ntt(a, b)
    if printPoly: print("c (cyclic NTT)=", c)
    print(f"equal: {c == c_red}")
    if c == c_red:
        return 0
    else:
        return 1

def testcase_negacyclic_ntt(n, q, printPoly=True):
    print(f"Testing naive negacyclic NTT multiplication (i.e., mod x^n + 1) with q={q}, n={n}")
    a = Poly.random(n, q)
    b = Poly.random(n, q)
    if printPoly: print("a=", a)
    if printPoly: print("b=", b)

    # compute reference product using schoolbook for comparison
    c_ref = a*b
    # reduce mod x^n + 1
    c_red = Poly(c_ref.coeffs[:n], q)
    for i in range(n, c_ref.n):
        c_red.coeffs[i-n] -= c_ref.coeffs[i]
    c_red.reduce()

    if printPoly: print("a*b (ref)=", c_red)
    c = polymul_negacyclic_ntt(a, b)
    if printPoly: print("c (cyclic NTT)=", c)
    print(f"equal: {c == c_red}")
    if c == c_red:
        return 0
    else:
        return 1

if __name__ == "__main__":
    rc = 0

    # plain NTT (no reduction)
    rc |= testcase_ntt(n=4, q=29)
    rc |= testcase_ntt(n=256, q=3067, printPoly=False)

    # cyclic NTT (mod x^n-1)
    rc |= testcase_cyclic_ntt(n=8, q=17)
    rc |= testcase_cyclic_ntt(n=256, q=3329, printPoly=False)

    # negacyclic NTT (mod x^n+1)
    rc |= testcase_negacyclic_ntt(n=8, q=17)
    rc |= testcase_negacyclic_ntt(n=256, q=7681, printPoly=False)

    if rc != 0:
        print("TEST FAILED.")
        sys.exit(1)
    print("ALL GOOD.")