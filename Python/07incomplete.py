#!/usr/bin/env python3
from poly import Poly
from common import bitreverse, isPrime, primitiveRootOfUnity
import sys
import math

# TODO: write documentation for each function
def precomp_ct_cyclic(n, root, q, numLayers):
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    assert numLayers <= logn

    twiddles = [pow(root, i, q) for i in range(2**numLayers//2)]
    twiddles = bitreverse(twiddles)

    twiddlesPerLayer = []
    for i in range(numLayers):
        twiddlesPerLayer.append(twiddles[:2**i])

    return twiddlesPerLayer

def precomp_basemul_cyclic(n, root, q, numLayers):
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    assert numLayers <= logn

    twiddles = [pow(root, i, q) for i in range(2**numLayers)]
    twiddles = bitreverse(twiddles)
    return twiddles

def precomp_gs_cyclic(n, root, q, numLayers):
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    twiddles = [pow(root, -i, q) for i in range(2**numLayers//2)]
    twiddles = bitreverse(twiddles)

    twiddlesPerLayer = []
    for i in range(logn-numLayers, logn):
        twiddlesPerLayer.append(twiddles[:2**(logn - 1- i)])

    return twiddlesPerLayer


def precomp_ct_negacyclic(n, root, q, numLayers):
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    twiddles = [pow(root, i, q) for i in range(2**numLayers)]
    twiddles = bitreverse(twiddles)
    twiddlesPerLayer = []
    off = 1
    for i in range(numLayers):
        twiddlesPerLayer.append(twiddles[off:off+2**i])
        off = off+2**i

    return twiddlesPerLayer

def precomp_basemul_negacyclic(n, root, q, numLayers):
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    assert numLayers <= logn

    twiddles = [pow(root, 2*i+1, q) for i in range(2**numLayers)]
    twiddles = bitreverse(twiddles)
    return twiddles


def precomp_gs_negacyclic(n, root, q, numLayers):
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    twiddles = [pow(root, -(i+1), q) for i in range(2**numLayers)]
    twiddles = bitreverse(twiddles)

    twiddlesPerLayer = []
    off = 0
    for i in range(logn-numLayers, logn):
        twiddlesPerLayer.append(twiddles[off:off+2**(logn-1-i)])
        off = off+2**(logn-1-i)

    return twiddlesPerLayer

def ntt_ct(a, twiddles, numLayers):
    a = a.copy()
    n = a.n
    q = a.q
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    assert len(twiddles) == numLayers

    def ct_butterfly(a0, a1, twiddle, q):
        tmp = (a1 * twiddle) % q
        return (a0 + tmp) % q, (a0 + q - tmp) % q

    for i in range(numLayers):
        distance = 2**(logn - 1 - i)
        for j in range(2**i):
            twiddle = twiddles[i][j]
            for k in range(distance):
                idx0 = 2*j*distance + k
                idx1 = idx0 + distance
                a.coeffs[idx0], a.coeffs[idx1] = ct_butterfly(a.coeffs[idx0], a.coeffs[idx1], twiddle, q)

    return a

def invntt_gs(a, twiddles, numLayers):
    a = a.copy()
    n = a.n
    q = a.q
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    assert len(twiddles) == numLayers

    def gs_butterfly(a0, a1, twiddle, q):
        tmp = (a0 + a1) % q
        a1  = ((a0 - a1)*twiddle) % q
        a0  = tmp
        return (a0, a1)

    for i in range(logn-numLayers, logn):
        distance = 2**(i)
        for j in range(2**(logn - 1 - i)):
            twiddle = twiddles[i-logn+numLayers][j]
            for k in range(distance):
                idx0 = 2*j*distance + k
                idx1 = idx0 + distance
                a.coeffs[idx0], a.coeffs[idx1] = gs_butterfly(a.coeffs[idx0], a.coeffs[idx1], twiddle, q)

    # Note that half of these multiplications can be merged in the last round
    # of GS butterflies, by precomputing (twiddle_i*n^-1) % q.
    ninv = pow(2**numLayers, -1, q)
    for i in range(a.n):
        a.coeffs[i] = (a.coeffs[i] * ninv) % q
    return a



def polymul_ntt_ct_gs_incomplete(a, b, twiddleNtt, twiddlesInvntt, twiddlesBasemul, numLayers):
    n = a.n
    q = a.q
    logn =  int(math.log(n, 2))
    pointwiseDegree = 2**(logn-numLayers)

    antt = ntt_ct(a, twiddleNtt, numLayers)
    bntt = ntt_ct(b, twiddleNtt, numLayers)

    cntt = Poly.zero(n, q)
    for i in range(a.n//pointwiseDegree):
        a = Poly(antt.coeffs[i*pointwiseDegree:(i+1)*pointwiseDegree], q)
        b = Poly(bntt.coeffs[i*pointwiseDegree:(i+1)*pointwiseDegree], q)
        c = a*b
        for j in range(pointwiseDegree):
            cntt.coeffs[i*pointwiseDegree+j] = c.coeffs[j]
        for j in range(pointwiseDegree, c.n):
            cntt.coeffs[i*pointwiseDegree+j-pointwiseDegree]  += c.coeffs[j]*twiddlesBasemul[i]

    cntt.reduce()
    return invntt_gs(cntt, twiddlesInvntt, numLayers)


def testcase_cyclic(n, q, numLayers, printPoly=True):
    rc = 0
    print(f"Testing polynomial multiplication using cyclic incomplete NTT ({numLayers} layers) with n={n}, q={q}")
    # find a 2**numLayers-th root of unity
    root = primitiveRootOfUnity(2**numLayers, q)

    # precompute twiddles
    twiddlesNtt = precomp_ct_cyclic(n, root, q, numLayers)
    twiddlesInvntt = precomp_gs_cyclic(n, root, q, numLayers)
    twiddlesBasemul = precomp_basemul_cyclic(n, root, q, numLayers)

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
    c = polymul_ntt_ct_gs_incomplete(a, b, twiddlesNtt, twiddlesInvntt, twiddlesBasemul, numLayers)
    if printPoly: print("c (cyclic NTT)=", c)
    print(f"equal: {c == c_red}")
    if c != c_red:
        rc = 1
    return rc

def testcase_negacyclic(n, q, numLayers, printPoly=True):
    rc = 0
    print(f"Testing polynomial multiplication using negacyclic incomplete NTT ({numLayers} layers) with n={n}, q={q}")
    # find a 2**(numLayers+1)-th root of unity
    root = primitiveRootOfUnity(2**(numLayers+1), q)

    # precompute twiddles
    twiddlesNtt = precomp_ct_negacyclic(n, root, q, numLayers)
    twiddlesInvntt = precomp_gs_negacyclic(n, root, q, numLayers)
    twiddlesBasemul = precomp_basemul_negacyclic(n, root, q, numLayers)

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
    c = polymul_ntt_ct_gs_incomplete(a, b, twiddlesNtt, twiddlesInvntt, twiddlesBasemul, numLayers)
    if printPoly: print("c (cyclic NTT)=", c)
    print(f"equal: {c == c_red}")
    if c != c_red:
        rc = 1
    return rc

if __name__ == "__main__":
    rc = 0

    # test cyclic
    rc |= testcase_cyclic(n=8, q=5, numLayers=2)
    rc |= testcase_cyclic(n=8, q=5, numLayers=1)
    rc |= testcase_cyclic(n=256, q=3329, numLayers=7, printPoly=False)
    rc |= testcase_cyclic(n=256, q=3329, numLayers=6, printPoly=False)

    # test negacyclic
    rc |= testcase_negacyclic(n=8, q=17, numLayers=2)
    rc |= testcase_negacyclic(n=8, q=17, numLayers=1)
    # The Kyber NTT
    rc |= testcase_negacyclic(n=256, q=3329, numLayers=7, printPoly=False)
    rc |= testcase_negacyclic(n=256, q=3329, numLayers=6, printPoly=False)

    if rc != 0:
        print("TEST FAILED.")
        sys.exit(1)
    print("ALL GOOD.")