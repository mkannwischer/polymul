#!/usr/bin/env python3
import sys
import math
from poly import Poly
# TODO: change all imports to not use specific imports; so that one know where to look for stuff
from common import isPrime, primitiveRootOfUnity, reverseBaseN, \
                   ntt_naive_cyclic, invntt_naive_cyclic

# TODO: write documentation for each function
def precomp_ntt_cyclic(n, root, q):
    logn =  int(math.log(n, 3))
    assert 3**logn == n

    twiddles = []
    for j in reverseBaseN(list(range(3**(logn-1))), 3):
        twiddles.extend([pow(root, j, q), pow(root, 2*j, q)])

    twiddlesPerLayer = []
    for i in range(logn):
        twiddlesPerLayer.append(twiddles[:2*3**i])

    print(twiddlesPerLayer)
    return twiddlesPerLayer

def precomp_invntt_cyclic(n, root, q):

    logn =  int(math.log(n, 3))
    assert 3**logn == n

    twiddles = []
    for j in reverseBaseN(list(range(3**(logn-1))), 3):
        twiddles.extend([pow(root, -j, q), pow(root, -2*j, q)])

    twiddlesPerLayer = []
    for i in range(logn):
        twiddlesPerLayer.append(twiddles[:2*3**(logn - 1- i)])

    print(twiddlesPerLayer)
    return twiddlesPerLayer

def ntt(a, twiddles, root3):
    a = a.copy()
    n = a.n
    q = a.q
    logn =  int(math.log(n, 3))
    assert 3**logn == n

    def butterfly(a0, a1, a2, twiddles, q, root3):
        # a_i     = a_i +     c a_{i+n} +     c^2 a_{i+2n}
        # a_{i+n} = a_i +   w c a_{i+n} + w^2 c^2 a_{i+2n}
        # a_{i+2n}= a_i + w^2 c a_{i+n} + w   c^2 a_{i+2n}
        # twiddles are [c, c^2]
        tmp0 = a0 + a1*twiddles[0]             + a2*twiddles[1]
        tmp1 = a0 + a1*twiddles[0]*root3       + a2*twiddles[1]*root3*root3
        tmp2 = a0 + a1*twiddles[0]*root3*root3 + a2*twiddles[1]*root3
        return (tmp0 % q, tmp1 % q, tmp2 % q)


    for i in range(logn):
        distance = 3**(logn - 1 - i)
        for j in range(3**i):
            # one butterfly requires 2 twiddles
            twiddleList = twiddles[i][j*2:(j+1)*2]
            for k in range(distance):
                idx0 = 3*j*distance + k
                idx1 = idx0 + distance
                idx2 = idx1 + distance
                a.coeffs[idx0], a.coeffs[idx1], a.coeffs[idx2] = butterfly(a.coeffs[idx0], a.coeffs[idx1], a.coeffs[idx2], twiddleList, q, root3)
    return a

def invntt(a, twiddles, root3):
    a = a.copy()
    n = a.n
    q = a.q
    logn =  int(math.log(n, 3))
    assert 3**logn == n

    def butterfly(a0, a1, a2, twiddles, q, root3):
        # a_i     = 1/3     (a_i +     a_{i+n} +     a_{i+2n})
        # a_{i+n} = 1/(3c)  (a_i + w^2 a_{i+n} + w   a_{i+2n})
        # a_{i+2n}= 1/(3c^2)(a_i + w   a_{i_b} + w^2 a_{i+2n})

        # We delay the multiplyications by 1/3 til the end
        tmp0 =             (a0 + a1 + a2)
        tmp1 = twiddles[0]*(a0 + a1*root3*root3 + a2*root3)
        tmp2 = twiddles[1]*(a0 + a1*root3       + a2*root3*root3)
        return (tmp0 % q, tmp1 % q, tmp2 % q)


    for i in range(logn):
        distance = 3**(i)
        for j in range(3**(logn - 1 - i)):
            # one butterfly requires 2 twiddles
            twiddleList = twiddles[i][j*2:(j+1)*2]
            for k in range(distance):
                idx0 = 3*j*distance + k
                idx1 = idx0 + distance
                idx2 = idx1 + distance
                a.coeffs[idx0], a.coeffs[idx1], a.coeffs[idx2] = butterfly(a.coeffs[idx0], a.coeffs[idx1], a.coeffs[idx2], twiddleList, q, root3)
    # divide by n^-1
    ninv = pow(a.n, -1, q)
    for i in range(a.n):
        a.coeffs[i] = (a.coeffs[i] * ninv) % q
    return a

def polymul_ntt(a, b, twiddleNtt, twiddlesInvntt, root3):
    antt = ntt(a, twiddleNtt, root3)
    bntt = ntt(b, twiddleNtt, root3)
    cntt = antt.pointwise(bntt)
    return invntt(cntt, twiddlesInvntt, root3)

def testcase_cyclic(n, q, printPoly=True):
    rc = 0
    logn = int(math.log(n, 3))
    # find an n-th root of unity and the corresponding 3-rd root of unity
    root  = primitiveRootOfUnity(n, q)
    root3 = pow(root, 3**(logn-1), q)

    # precompute twiddles
    twiddlesNtt = precomp_ntt_cyclic(n, root, q)
    twiddlesInvNtt = precomp_invntt_cyclic(n, root, q)

    print(f"Testing forward cyclic NTT using CT butterflies with n={n}, q={q}")
    # pick a random polynomial for testing
    a = Poly.random(n, q)
    if printPoly: print("a=", a)
    # tranform to ntt domain (output will be base-3-reversed)
    antt = ntt(a, twiddlesNtt, root3)
    if printPoly: print("antt=", antt)
    # compute a naive reference NTT to compare (output is in normal order)
    anttref = ntt_naive_cyclic(a, root)
    # base-3-reverse, so we can compare
    anttref.coeffs = reverseBaseN(anttref.coeffs, 3)
    if printPoly: print("antt_ref=", anttref)
    print(f"equal: {antt == anttref}")
    if antt != anttref:
        rc = 1

    print(f"Testing inverse cyclic NTT using GS butterflies with n={n}, q={q}")
    # pick a random polynomial for testing
    antt = Poly.random(n, q)
    if printPoly: print("antt=", antt)
    # GS needs inputs in base-3-reversed order
    anttbrv = Poly(reverseBaseN(antt.coeffs, 3), q)
    # tranform to ntt domain (output will be base-3-reversed)
    a = invntt(anttbrv, twiddlesInvNtt, root3)
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
    antt = ntt(a, twiddlesNtt, root3)
    if printPoly: print("antt=", antt)
    # compute inverse transform
    a2  = invntt(antt, twiddlesInvNtt, root3)
    if printPoly: print("a2=", a2)
    # check that a ==invntt(ntt(a))
    print(f"equal: {a == a2}")
    if a != a2:
        rc = 1

    print(f"Testing polynomial multiplication cyclic with n={n}, q={q}")
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
    c = polymul_ntt(a, b, twiddlesNtt, twiddlesInvNtt, root3)
    if printPoly: print("c (cyclic NTT)=", c)
    print(f"equal: {c == c_red}")
    if c != c_red:
        rc = 1
    return rc


if __name__ == "__main__":
    rc = 0

    rc |= testcase_cyclic(n=3, q=7)
    rc |= testcase_cyclic(n=9, q=37)
    rc |= testcase_cyclic(n=27, q=109)

    if rc != 0:
        print("TEST FAILED.")
        sys.exit(1)
    print("ALL GOOD.")