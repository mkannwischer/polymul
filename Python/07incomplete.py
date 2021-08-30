#!/usr/bin/env python3
"""Section 2.2.7: Incomplete NTT.

Covers fast fourier-transforms implementing incomplete NTTs.
"""
from poly import Poly
import common
import sys
import math


def precomp_ct_cyclic(n, root, q, numLayers):
    """Precompute the required twiddle factors for a incomplete cyclic Cooley--Tukey FFT.

    First layer: [1]
    Second layer: [1, -1] = [1, root^(n/2)]
    Third layer: [1, -1, sqrt(-1), -sqrt(-1)] = [1, root^(n/2), root^(n/4), root^(3n/4)]
    ...

    Note that the twiddle factors repeat. In a real implementation one would
    not store them repeatedly. One would also eliminate the multiplications
    by 1 and -1.

    Parameters
    ----------
    n : int
        size of the NTT (number of coefficients).
    root : int
        (2^numLayers)-th primitive root of unity modulo q.
    q : int
        modulus.
    numLayers: int
        number of layers in the NTT. Needs to be <= log n.
    Returns
    ----------
    list
        list of lists of twiddle factors.
    """
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    assert numLayers <= logn

    twiddles = [pow(root, i, q) for i in range(2**numLayers//2)]
    twiddles = common.bitreverse(twiddles)

    twiddlesPerLayer = []
    for i in range(numLayers):
        twiddlesPerLayer.append(twiddles[:2**i])

    return twiddlesPerLayer

def precomp_basemul_cyclic(n, root, q, numLayers):
    """Precompute the required twiddle factors for the base multiplication of an incomplete cyclic NTT.

    Parameters
    ----------
    n : int
        size of the NTT (number of coefficients).
    root : int
        (2^numLayers)-th primitive root of unity modulo q.
    q : int
        modulus.
    numLayers: int
        number of layers in the NTT. Needs to be <= log n.
    Returns
    ----------
    list
        list of lists of twiddle factors.
    """
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    assert numLayers <= logn

    twiddles = [pow(root, i, q) for i in range(2**numLayers)]
    twiddles = common.bitreverse(twiddles)
    return twiddles

def precomp_gs_cyclic(n, root, q, numLayers):
    """Precompute the required twiddle factors for a incomplete cyclic Gentleman--Sande inverse FFT.

    The twiddles correspond to the inverses of the ones computes in `precomp_ct_cyclic`.
    Note that the twiddle factors repeat. In a real implementation one would
    not store them repeatedly.

    Parameters
    ----------
    n : int
        size of the NTT (number of coefficients).
    root : int
        (2^numLayers)-th primitive root of unity modulo q.
    q : int
        modulus.
    numLayers: int
        number of layers in the NTT. Needs to be <= log n.
    Returns
    ----------
    list
        list of lists of twiddle factors.
    """
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    twiddles = [pow(root, -i, q) for i in range(2**numLayers//2)]
    twiddles = common.bitreverse(twiddles)

    twiddlesPerLayer = []
    for i in range(logn-numLayers, logn):
        twiddlesPerLayer.append(twiddles[:2**(logn - 1- i)])

    return twiddlesPerLayer


def precomp_ct_negacyclic(n, root, q, numLayers):
    """Precompute the required twiddle factors for a incomplete negacyclic Cooley--Tukey FFT.

    First layer: [-1] = [root^(n/2)]
    Second layer: [sqrt(-1), -sqrt(-1)] = [root^(n/4), root^(3n/4)]
    Third layer: [sqrt(root^(n/4)), -sqrt(root^(n/4)), sqrt(root^(3n/4)), -sqrt(root^(3n/4))]
                =[root^(n/8), root^(5n/8), root^(3n/8), root^(7n/8)]
    ...

    Parameters
    ----------
    n : int
        size of the NTT (number of coefficients).
    root : int
        2*(2^numLayers)-th primitive root of unity modulo q.
    q : int
        modulus.
    numLayers: int
        number of layers in the NTT. Needs to be <= log n.
    Returns
    ----------
    list
        list of lists of twiddle factors.
    """
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    twiddles = [pow(root, i, q) for i in range(2**numLayers)]
    twiddles = common.bitreverse(twiddles)
    twiddlesPerLayer = []
    off = 1
    for i in range(numLayers):
        twiddlesPerLayer.append(twiddles[off:off+2**i])
        off = off+2**i

    return twiddlesPerLayer

def precomp_basemul_negacyclic(n, root, q, numLayers):
    """Precompute the required twiddle factors for the base multiplication of an incomplete negacyclic NTT.

    Parameters
    ----------
    n : int
        size of the NTT (number of coefficients).
    root : int
        2*(2^numLayers)-th primitive root of unity modulo q.
    q : int
        modulus.
    numLayers: int
        number of layers in the NTT. Needs to be <= log n.
    Returns
    ----------
    list
        list of lists of twiddle factors.
    """
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    assert numLayers <= logn

    twiddles = [pow(root, 2*i+1, q) for i in range(2**numLayers)]
    twiddles = common.bitreverse(twiddles)
    return twiddles


def precomp_gs_negacyclic(n, root, q, numLayers):
    """Precompute the required twiddle factors for a incomplete negacyclic Gentleman--Sande inverse FFT.

    The twiddles correspond to the inverses of the ones computes in `precomp_ct_negacyclic`.

    Parameters
    ----------
    n : int
        size of the NTT (number of coefficients).
    root : int
        2*(2^numLayers)-th primitive root of unity modulo q.
    q : int
        modulus.
    numLayers: int
        number of layers in the NTT. Needs to be <= log n.
    Returns
    ----------
    list
        list of lists of twiddle factors.
    """
    logn =  int(math.log(n, 2))
    assert 2**logn == n
    twiddles = [pow(root, -(i+1), q) for i in range(2**numLayers)]
    twiddles = common.bitreverse(twiddles)

    twiddlesPerLayer = []
    off = 0
    for i in range(logn-numLayers, logn):
        twiddlesPerLayer.append(twiddles[off:off+2**(logn-1-i)])
        off = off+2**(logn-1-i)

    return twiddlesPerLayer

def ntt_ct(a, twiddles, numLayers):
    """Compute a Cooley--Tukey FFT. Stop after numLayers.

    Expects twiddles to be computed by `precomp_ct_cyclic` or `precomp_ct_negacyclic`
    Each layer computes a split of
    Z_q[x]/(x^n - c^2) to Z_q[x]/(x^(n/2) - c) x Z_q[x]/(x^(n/2) + c)
    using the CT butterfly:
        a_i' = a_i + c*a_j
        a_j' = a_i - c*a_j

    Parameters
    ----------
    a : list
        polynomial with n coefficients to be transformed to NTT domain.
    twiddles : list
        list of lists of twiddle factors per NTT layer.
    numLayers: list
        number of layers in the NTT. Needs to be <= log n.
    Returns
    ----------
    Poly
        a transformed to NTT domain.
    """
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
    """Compute a Gentleman--Sande inverse FFT. Stop after numLayers.

    Expects twiddles to be computed by `precomp_gs_cyclic` or `precomp_gs_negacyclic`
    Each layer computes the CRT of
    Z_q[x]/(x^(n/2) - c) x Z_q[x]/(x^(n/2) + c) to recover an element in Z_q[x]/(x^n - c^2)
    using the GS butterfly:
        a_i' = 1/2 * (a_i + a_j)
        a_j' = 1/2 * 1/c * (a_i - a_j)
    The scaling by 1/2 is usually delayed until the very end, i.e., multiplication by 1/(2^numLayers).

    Parameters
    ----------
    a : list
        input in NTT domain.
    twiddles : list
        list of lists of twiddle factors per NTT layer.
    numLayers: list
        number of layers in the NTT. Needs to be <= log n.
    Returns
    ----------
    Poly
        a transformed to normal domain.
    """
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
    """Compute a polynomial multiplication by computing iNTT(NTT(a) o NTT(b)) using incomplete NTTs and o denoting basemul.

    Works for both the cyclic and the negacyclic case (with the correct twiddles).

    Parameters
    ----------
    a : Poly
        first multiplicand polynomial with n coefficients.
    b : Poly
        second multiplicand polynomial with n coefficients.
    twiddlesNtt : list
        twiddles for the foward NTT as computed by `precomp_ct_cyclic` or `precomp_ct_negacyclic`.
    tiwddlesInvntt : list
        twiddles for the inverse nTT as computed by `precomp_gs_cyclic` or `precmp_gs_negacyclic`.
    numLayers: list
        number of layers in the NTT. Needs to be <= log n.
    Returns
    ----------
    Poly
        product a*b with n coefficients.
    """
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
    """Random test of cyclic NTT multiplication for Zq[x]/(x^n-1).

    Parameters
    ----------
    n : int
        number of coefficients of input polynomials.
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
    print(f"Testing polynomial multiplication using cyclic incomplete NTT ({numLayers} layers) with n={n}, q={q}")
    # find a 2**numLayers-th root of unity
    root = common.primitiveRootOfUnity(2**numLayers, q)

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
    """Random test of negacyclic NTT multiplication for Zq[x]/(x^n+1).

    Parameters
    ----------
    n : int
        number of coefficients of input polynomials.
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
    print(f"Testing polynomial multiplication using negacyclic incomplete NTT ({numLayers} layers) with n={n}, q={q}")
    # find a 2**(numLayers+1)-th root of unity
    root = common.primitiveRootOfUnity(2**(numLayers+1), q)

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