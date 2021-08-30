#!/usr/bin/env python3
"""Section 2.2.1: Schoolbook multiplication.

Illustrates schoolbook multiplication in various polynomial rings.
This is not rocket science, but for reference of what the other algorithms
should be implementing.

"""

from poly import Poly


def polymul_schoolbook(a, b):
    """Multiplies two polynomials a and b in Z_q[x].

    Parameters
    ----------
    a : Poly
        first multiplicand with n1 coefficients.
    b : Poly
        second multiplicand with n2 coefficients.
    Returns
    ----------
    Poly
        product a*b with n1+n2-1 coefficients.
    """
    assert a.q == b.q

    # if multiplicands are of degree at most  n, the product will be of degree at most n+n
    # i.e., if multiplicands have n+1 coefficients, the product will have n+n+1 coefficients
    c = Poly.zero(a.n+b.n-1, a.q)
    for i in range(a.n):
        for j in range(b.n):
            c.coeffs[i+j] += a.coeffs[i]*b.coeffs[j]
    # reduce coefficients modulo q
    if c.q != None:
        for i in range(c.n):
            c.coeffs[i] %= c.q
    return c

def polymul_schoolbook_cyclic(a, b):
    """Cyclic convolution of two polynomials, i.e., multiplication mod x^n - 1.

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
    assert a.q == b.q
    assert a.n == b.n
    n = a.n
    q = a.q
    # compute product with 2*n-1 coefficients first
    c2n = a*b
    # then perform convolution x^N = 1
    c = Poly(c2n.coeffs[:n], q)
    for i in range(n, c2n.n):
        c.coeffs[i-n] += c2n.coeffs[i]

    c.reduce()
    return c

def polymul_schoolbook_negacyclic(a, b):
    """Negacyclic convolution of two polynomials, i.e., multiplication mod x^n + 1.

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
    assert a.q == b.q
    assert a.n == b.n
    n = a.n
    q = a.q
    # compute product with 2*n-1 coefficients first
    c2n = a*b
    # then perform convolution x^N = -1
    c = Poly(c2n.coeffs[:n], q)
    for i in range(n, c2n.n):
        c.coeffs[i-n] -= c2n.coeffs[i]

    c.reduce()
    return c

if __name__ == "__main__":
    a = Poly([0, 4])
    b = Poly([13, 4])
    print("a=", a) # a= 4x^1
    print("b=", b) # b= 4x^1 + 13x^0

    print("# Polynomial multiplication in Z[x]")
    c = polymul_schoolbook(a, b)
    print("a*b=", c) # a*b= 16x^2 + 52x^1

    print("# Polynomial multiplication in Z_{17}[x]")
    a = Poly([0, 4], q=17)
    b = Poly([13, 4], q=17)
    c = polymul_schoolbook(a, b)
    print("a*b=", c) # a*b= 16x^2 + 1x^1

    # for convienence, the Poly class can also do that for us:
    print("a*b=", a*b) # a*b= 16x^2 + 1x^1

    print("# Polynomial multiplication in Z_{17}[x]/(x^2-1)")
    c = polymul_schoolbook_cyclic(a, b)
    print("a*b=", c) # a*b= 1x^1 + 16x^0

    print("# Polynomial multiplication in Z_{17}[x]/(x^2+1)")
    c = polymul_schoolbook_negacyclic(a, b)
    print("a*b=", c) # a*b= 1x^1 + 1x^0
