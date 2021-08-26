#!/usr/bin/env python3
from poly import Poly
from common import bitreverse, isPrime, primitiveRootOfUnity, ntt_naive_cyclic, \
                   invntt_naive_cyclic, ntt_naive_negacyclic, invntt_naive_negacyclic
import sys
import math

# TODO: write documentation for each function

def precomp_ct_cyclic(n, root, q):
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    twiddles = [pow(root, i, q) for i in range(n//2)]
    twiddles = bitreverse(twiddles)

    twiddlesPerLayer = []
    for i in range(logn):
        twiddlesPerLayer.append(twiddles[:2**i])

    return twiddlesPerLayer


def precomp_ct_negacyclic(n, root, q):
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    twiddles = [pow(root, i, q) for i in range(n)]
    twiddles = bitreverse(twiddles)
    twiddlesPerLayer = []
    off = 1
    for i in range(logn):
        twiddlesPerLayer.append(twiddles[off:off+2**i])
        off = off+2**i

    return twiddlesPerLayer

def ntt_ct(a, twiddles):
    a = a.copy()
    n = a.n
    q = a.q
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    assert len(twiddles) == logn

    def ct_butterfly(a0, a1, twiddle, q):
        tmp = (a1 * twiddle) % q
        return (a0 + tmp) % q, (a0 + q - tmp) % q

    for i in range(logn):
        distance = 2**(logn - 1 - i)
        for j in range(2**i):
            twiddle = twiddles[i][j]
            for k in range(distance):
                idx0 = 2*j*distance + k
                idx1 = idx0 + distance
                a.coeffs[idx0], a.coeffs[idx1] = ct_butterfly(a.coeffs[idx0], a.coeffs[idx1], twiddle, q)

    return a

def precomp_gs_cyclic(n, root, q):
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    twiddles = [pow(root, -i, q) for i in range(n//2)]
    twiddles = bitreverse(twiddles)

    twiddlesPerLayer = []
    for i in range(logn):
        twiddlesPerLayer.append(twiddles[:2**(logn - 1- i)])

    return twiddlesPerLayer

def precomp_gs_negacyclic(n, root, q):
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    twiddles = [pow(root, -(i+1), q) for i in range(n)]
    twiddles = bitreverse(twiddles)

    twiddlesPerLayer = []
    off = 0
    for i in range(logn):
        twiddlesPerLayer.append(twiddles[off:off+2**(logn-1-i)])
        off = off+2**(logn-1-i)

    return twiddlesPerLayer

def invntt_gs(a, twiddles):
    a = a.copy()
    n = a.n
    q = a.q
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    assert len(twiddles) == logn

    def gs_butterfly(a0, a1, twiddle, q):
        tmp = (a0 + a1) % q
        a1  = ((a0 - a1)*twiddle) % q
        a0  = tmp
        return (a0, a1)

    for i in range(logn):
        distance = 2**(i)
        for j in range(2**(logn - 1 - i)):
            twiddle = twiddles[i][j]
            for k in range(distance):
                idx0 = 2*j*distance + k
                idx1 = idx0 + distance
                a.coeffs[idx0], a.coeffs[idx1] = gs_butterfly(a.coeffs[idx0], a.coeffs[idx1], twiddle, q)

    # Note: In real implementations half of these multiplications can be merged
    # into the last layer of the GS butterflies.
    ninv = pow(a.n, -1, q)
    for i in range(a.n):
        a.coeffs[i] = (a.coeffs[i] * ninv) % q
    return a

def polymul_ntt_ct_gs(a, b, twiddleNtt, twiddlesInvntt):
    antt = ntt_ct(a, twiddleNtt)
    bntt = ntt_ct(b, twiddleNtt)
    cntt = antt.pointwise(bntt)
    return invntt_gs(cntt, twiddlesInvntt)



def testcase_cyclic(n, q, printPoly=True):
    rc = 0
    # find a n-th root of unity
    root = primitiveRootOfUnity(n, q)
    # precompute twiddles
    twiddlesNtt = precomp_ct_cyclic(n, root, q)
    # precompute twiddles
    twiddlesInvntt = precomp_gs_cyclic(n, root, q)

    print(f"Testing forward cyclic NTT using CT butterflies with n={n}, q={q}")
    # pick a random polynomial for testing
    a = Poly.random(n, q)
    if printPoly: print("a=", a)
    # tranform to ntt domain (output will be bitreversed)
    antt = ntt_ct(a, twiddlesNtt)
    if printPoly: print("antt=", antt)
    # compute a naive reference NTT to compare (output is in normal order)
    anttref = ntt_naive_cyclic(a, root)
    # bitreverse, so we can compare
    anttref.coeffs = bitreverse(anttref.coeffs)
    if printPoly: print("antt_ref=", anttref)
    print(f"equal: {antt == anttref}")
    if antt != anttref:
        rc = 1

    print(f"Testing inverse cyclic NTT using GS butterflies with n={n}, q={q}")
    # pick a random polynomial for testing
    antt = Poly.random(n, q)
    if printPoly: print("antt=", antt)
    # GS needs inputs in bitreversed order
    anttbrv = Poly(bitreverse(antt.coeffs), q)
    # tranform to ntt domain (output will be bitreversed)
    a = invntt_gs(anttbrv, twiddlesInvntt)
    if printPoly: print("a=", a)
    # compute a naive reference NTT to compare (output is in normal order)
    aref = invntt_naive_cyclic(antt, root)
    if printPoly: print("aref=", aref)
    print(f"equal: {a == aref}")
    if a != aref:
        rc = 1

    print(f"Testing inverse a == invntt(ntt(a)) cyclic with n={n}, q={q}")
    a = Poly.random(n, q)
    if printPoly: print("a=", a)
    # compute forward transform
    antt = ntt_ct(a, twiddlesNtt)
    if printPoly: print("antt=", antt)
    # compute inverse transform
    a2  = invntt_gs(antt, twiddlesInvntt)
    if printPoly: print("a2=", a2)
    # check that a ==invntt(ntt(a))
    print(f"equal: {a == a2}")
    if a != a2:
        rc = 1

    print(f"Testing polynomial multiplication cyclic with n={n}, q={q}")
    a = Poly.random(n, q)
    b = Poly.random(n, q)
    print("a=", a)
    print("b=", b)

    # compute reference product using schoolbook for comparison
    c_ref = a*b
    # reduce mod x^n - 1
    c_red = Poly(c_ref.coeffs[:n], q)
    for i in range(n, c_ref.n):
        c_red.coeffs[i-n] += c_ref.coeffs[i]
    c_red.reduce()

    if printPoly: print("a*b (ref)=", c_red)
    c = polymul_ntt_ct_gs(a, b, twiddlesNtt, twiddlesInvntt)
    if printPoly: print("c (cyclic NTT)=", c)
    print(f"equal: {c == c_red}")
    if c != c_red:
        rc = 1
    return rc

def testcase_negacyclic(n, q, printPoly=True):
    rc = 0
    # find a n-th root of unity
    root = primitiveRootOfUnity(2*n, q)
    # precompute twiddles
    twiddlesNtt = precomp_ct_negacyclic(n, root, q)
    # precompute twiddles
    twiddlesInvntt = precomp_gs_negacyclic(n, root, q)

    print(f"Testing forward negacyclic NTT using CT butterflies with n={n}, q={q}")
    # pick a random polynomial for testing
    a = Poly.random(n, q)
    if printPoly: print("a=", a)
    # tranform to ntt domain (output will be bitreversed)
    antt = ntt_ct(a, twiddlesNtt)
    if printPoly: print("antt=", antt)
    # compute a naive reference NTT to compare (output is in normal order)
    anttref = ntt_naive_negacyclic(a, root)
    # bitreverse, so we can compare
    anttref.coeffs = bitreverse(anttref.coeffs)
    if printPoly: print("antt_ref=", anttref)
    print(f"equal: {antt == anttref}")
    if antt != anttref:
        rc = 1

    print(f"Testing inverse negacyclic NTT using GS butterflies with n={n}, q={q}")

    # pick a random polynomial for testing
    antt = Poly.random(n, q)
    if printPoly: print("antt=", antt)
    # GS needs inputs in bitreversed order
    anttbrv = Poly(bitreverse(antt.coeffs), q)
    # tranform to ntt domain (output will be bitreversed)
    a = invntt_gs(anttbrv, twiddlesInvntt)
    if printPoly: print("a=", a)
    # compute a naive reference NTT to compare (output is in normal order)
    aref = invntt_naive_negacyclic(antt, root)
    if printPoly: print("aref=", aref)
    print(f"equal: {a == aref}")
    if a != aref:
        rc = 1

    print(f"Testing inverse a == invntt(ntt(a)) negacyclic with n={n}, q={q}")
    a = Poly.random(n, q)
    if printPoly: print("a=", a)
    # compute forward transform
    antt = ntt_ct(a, twiddlesNtt)
    if printPoly: print("antt=", antt)
    # compute inverse transform
    a2  = invntt_gs(antt, twiddlesInvntt)
    if printPoly: print("a2=", a2)
    # check that a ==invntt(ntt(a))
    print(f"equal: {a == a2}")
    if a != a2:
        rc = 1

    print(f"Testing polynomial multiplication negacyclic with n={n}, q={q}")
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
    c = polymul_ntt_ct_gs(a, b, twiddlesNtt, twiddlesInvntt)
    if printPoly: print("c (cyclic NTT)=", c)
    print(f"equal: {c == c_red}")
    if c != c_red:
        rc = 1
    return rc

if __name__ == "__main__":
    rc = 0

    # test cyclic NTT (mod x^n-1)
    rc |= testcase_cyclic(n=8, q=17)
    rc |= testcase_cyclic(n=256, q=3329, printPoly=False)

    # test negacyclic NTT (mod x^n+1)
    rc |= testcase_negacyclic(n=8, q=17)
    rc |= testcase_negacyclic(n=256, q=7681, printPoly=False)

    if rc != 0:
        print("TEST FAILED.")
        sys.exit(1)
    print("ALL GOOD.")