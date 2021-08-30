#!/usr/bin/env python3
"""Section 2.2.4: Number-theoretic transform.

Illustrates the cyclic and negacyclic number-theoric transform.
This is a reference implementation that requires time O(n^2) for the transform.
For fast implementations, see `05fft.py`
"""
from poly import Poly
import common
import sys


def polymul_ntt(a, b):
    """Perform NTT-based multiplication of two polynomials in the ring Zq[x].

    For computing the product we evaluate the n-coefficient polynomials and b
    at (2n-1) powers of a (2n-1)-th root of unity, then multply coefficient wise,
    and interpolate the product using the inverse NTT.
    Primitive (2n-1)-th root of unity modulo q needs to exist.

    Parameters
    ----------
    a : Poly
        first multiplicand with n coefficients.
    b : Poly
        second multiplicand with n coefficients.
    Returns
    ----------
    Poly
        product a*b with 2n-1 coefficients.
    """
    assert a.n == b.n
    assert a.q == b.q
    n = a.n
    q = a.q
    assert common.isPrime(q)
    assert ((q-1) % (2*n-1)) == 0

    # For implementing a general purpose NTT (not cyclic, negacyclic), we need
    # to evaluate a and b at 2n-1 roots of unity.

    def ntt_naive(p, root):
        """Compute forward (2n-1) NTT of n-coefficient polynomial.

        Compute pntt_i = sum_{j=0}^{n} p[j]*root^{ij} for 0 <= i < 2n-1

        Parameters
        ----------
        p : Poly
            input polynomial with n coefficients.
        root: int
            (2n-1)-th root of unity modulo q.
        Returns
        ----------
        Poly
            p evaluated at root^0, ..., root^(2n-2), i.e., p transformed into NTT domain.
        """
        pntt = Poly.zero(2*p.n-1, q)
        for i in range(2*p.n-1):
            for j in range(p.n):
                pntt.coeffs[i] = (pntt.coeffs[i] +  p.coeffs[j]*pow(root, i*j, q)) % q
        return pntt

    def invntt_naive(pntt, root):
        """Compute inverse NTT.

        Note that the inputs here are in NTT domain, i.e., have (2n-1) coefficients.

        Compute p_i =  1/n * sum_{j=0}^{2n-1} pntt[j]*root^{-ij} for 0 <= i < 2n-1

        Parameters
        ----------
        pntt : Poly
            input in NTT domain, i.e., evaluation of p at root^0, ..., root^(2n-2).
        root: int
            (2n-1)-th root of unity modulo q.
        Returns
        ----------
        Poly
            p in normal domain.
        """
        # pntt is already padded, i.e., pntt.n is 2*a.n-1
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
    root = common.primitiveRootOfUnity(2*n-1, q)

    # transform a and b to NTT domain
    antt = ntt_naive(a, root)
    bntt = ntt_naive(b, root)

    # point-wise multiplication
    cntt = antt.pointwise(bntt)

    # trnasform the result back to normal domain
    c = invntt_naive(cntt, root)

    return c


def polymul_cyclic_ntt(a, b):
    """Perform NTT-based multiplication of two polynomials in the ring Zq[x]/(x^n-1), i.e., cyclic.

    For computing the product we evaluate the n-coefficient polynomials and b
    at n powers of a n-th root of unity, then multply coefficient wise,
    and interpolate the product using the inverse NTT.
    Primitive n-th root of unity modulo q needs to exist.

    Parameters
    ----------
    a : Poly
        first multiplicand with n coefficients.
    b : Poly
        second multiplicand with n coefficients.
    Returns
    ----------
    Poly
        product a*b with n coefficients.
    """
    assert a.n == b.n
    assert a.q == b.q
    n = a.n
    q = a.q
    assert common.isPrime(q)
    assert ((q-1) % (n)) == 0


    # we need a n-th root of unity
    root = common.primitiveRootOfUnity(n, q)

    # For implementing a cyclic NTT we need to evaluate a and b at n roots of unity.

    def ntt_naive(p, root):
        """Compute forward cyclic n-NTT of n-coefficient polynomial.

        Compute pntt_i = sum_{j=0}^{n} p[j]*root^{ij} for 0 <= i < n

        Parameters
        ----------
        p : Poly
            input polynomial with n coefficients.
        root: int
            n-th root of unity modulo q.
        Returns
        ----------
        Poly
            p evaluated at root^0, ..., root^n-1, i.e., p transformed into NTT domain.
        """
        q = p.q
        pntt = Poly.zero(p.n, q)
        for i in range(p.n):
            for j in range(p.n):
                pntt.coeffs[i] = (pntt.coeffs[i] +  p.coeffs[j]*pow(root, i*j, q)) % q
        return pntt

    def invntt_naive(pntt, root):
        """Compute cyclic inverse NTT.

        Compute p_i =  1/n * sum_{j=0}^{n} pntt[j]*root^{-ij} for 0 <= i < n

        Parameters
        ----------
        pntt : Poly
            input in NTT domain, i.e., evaluation of p at root^0, ..., root^n-1.
        root: int
            n-th root of unity modulo q.
        Returns
        ----------
        Poly
            p in normal domain.
        """
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
    """Perform NTT-based multiplication of two polynomials in the ring Zq[x]/(x^n+1), i.e., negacyclic.

    For computing the product, we first twist a and b to Zq[x]/(x^n-1) by
    multiplying by powers of the 2n-th root of unity,
    then evaluate at n powers of a n-th root of unity,
    then multply coefficient-wise,
    then interpolate the product using the inverse NTT,
    then twist back to Zq[x]/(x^n+1) by multiplying by power of the inverse
    2n-th root of unity.
    Primitive 2n-th (and consequently, n-th) root of unity modulo q needs to exist.

    Parameters
    ----------
    a : Poly
        first multiplicand with n coefficients.
    b : Poly
        second multiplicand with n coefficients.
    Returns
    ----------
    Poly
        product a*b with n coefficients.
    """
    assert a.n == b.n
    assert a.q == b.q
    n = a.n
    q = a.q
    assert common.isPrime(q)
    # We need a 2n-th root of unity rather than an n-th root of unity
    assert ((q-1) % (2*n)) == 0

    # we need a n-th root of unity
    root = common.primitiveRootOfUnity(2*n, q)

    def ntt_naive(p, root):
        """Compute forward negacyclic n-NTT of n-coefficient polynomial.

           Compute pntt_i = sum_{j=0}^{n} p[j]*root^{j}*root^{2ij} for 0 <= i < n.
        This is the same as first twisting to Zq[x]/(x^n-1) and then performing a cyclic NTT.

        Parameters
        ----------
        p : Poly
            input polynomial with n coefficients.
        root: int
            2n-th root of unity modulo q.
        Returns
        ----------
        Poly
            p in NTT domain.
        """
        q = p.q
        pntt = Poly.zero(p.n, q)
        for i in range(p.n):
            for j in range(p.n):
                pntt.coeffs[i] = (pntt.coeffs[i] +  p.coeffs[j]*pow(root, j, q)*pow(root, 2*i*j, q)) % q
        return pntt

    def invntt_naive(pntt, root):
        """Compute negacyclic inverse NTT.

        Compute p_i =  1/n * root^{-i} * sum_{j=0}^{n} pntt[j]*root^{-2ij} for 0 <= i < n

        Parameters
        ----------
        pntt : Poly
            input in NTT domain.
        root: int
            2n-th root of unity modulo q.
        Returns
        ----------
        Poly
            p in normal domain.
        """
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
    """Random test of NTT multiplication for Zq[x] (`polymul_ntt`).

    Parameters
    ----------
    n : int
        number of coefficients of input polynomials.
    q : int
        modulus.
    printPoly : boolean
        flag for printing inputs and outputs.
    Returns
    ----------
    int
        0 if test is successful, 1 otherwise.
    """
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
    """Random test of NTT multiplication for Zq[x]/(x^n-1) (`polymul_cyclic_ntt`).

    Parameters
    ----------
    n : int
        number of coefficients of input polynomials.
    q : int
        modulus.
    printPoly : boolean
        flag for printing inputs and outputs.
    Returns
    ----------
    int
        0 if test is successful, 1 otherwise.
    """
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
    """Random test of NTT multiplication for Zq[x]/(x^n+1) (`polymul_negacyclic_ntt`).

    Parameters
    ----------
    n : int
        number of coefficients of input polynomials.
    q : int
        modulus.
    printPoly : boolean
        flag for printing inputs and outputs.
    Returns
    ----------
    int
        0 if test is successful, 1 otherwise.
    """
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