#!/usr/bin/env python3
from poly import Poly
import sys

def polymul_karatsuba(a, b):
    """
        This computes a polynomial product using two-way Karatsuba.
        Assumes n is the same for both a and b.
        Works for odd n.

        Smaller polynomial multiplications are implemented using schoolbook.
        For example, for n=4:
        Let c = a0b0; d = (a0+a1)(b0+b1); e = a1b1;

        We compute:
          0    1    2    3    4    5    6
         ---- ---- ----      ---- ---- ----
        | c0 | c1 | c2 |    | e0 | e1 | e2 |
         ---- ---- ---- ---- ---- ----  ----
                + | d0 | d1 | d2 |
                   ---- ---- ----
                - | c0 | c1 | c2 |
                   ---- ---- ----
                - | e0 | e1 | e2 |
                   ---- ---- ----
    """
    assert a.n == b.n
    assert a.q == b.q
    n = a.n
    q = a.q

    nhalf = n//2

    a0 = Poly(a.coeffs[:nhalf], q)
    a1 = Poly(a.coeffs[nhalf:], q)
    b0 = Poly(b.coeffs[:nhalf], q)
    b1 = Poly(b.coeffs[nhalf:], q)

    # Y = x^{n//2}
    # a = a_1 + Y a_0
    # b = b_1 + Y b_1
    # a*b = (a_1 + Y a_0)*(b_0 + Y b_1)
    #     = Y^2 a_1 b_1 + Y (a_0 b_1 + a_1 b_0) + a_0 b_0
    #     = Y^2 a_1 b_1 + Y ((a_0 + a_1)(b_0 + b_1) - a_1 b_1 - a_0 b_0) + a_0 b_0

    a0b0 = a0*b0
    a1b1 = a1*b1

    # compute t = ((a_0 + a_1)(b_0 + b_1) - a_1 b_1 - a_0 b_0)
    a0a1 = a0+a1
    b0b1 = b0+b1
    t = a0a1*b0b1
    t -= a0b0
    t -= a1b1

    # compute a*b = Y^2 a_1 b_1 + Y t + a_0 b_0
    c = Poly.zero(2*n-1, q)
    for i in range(a0b0.n):
        c.coeffs[i] = a0b0.coeffs[i]

    for i in range(t.n):
        c.coeffs[i+nhalf] += t.coeffs[i]

    for i in range(a1b1.n):
        c.coeffs[i+2*nhalf] += a1b1.coeffs[i]

    # reduce coefficients mod q
    c.reduce()
    return c

def polymul_karatsuba_recursive(a, b, threshold):
    """
        Similar as polymul_karatsuba, but the smaller polynomial multiplications
        are also using karatsuba down to a threadhold from which we switch to
        schoolbook.
        Assumes n is the same for both a and b.
    """
    assert a.n == b.n
    assert a.q == b.q
    n = a.n
    q = a.q

    # if polynomials are small enough, just use schoolbook multiplicaton
    if a.n <= threshold:
        return a*b

    nhalf = n//2

    a0 = Poly(a.coeffs[:nhalf], q)
    a1 = Poly(a.coeffs[nhalf:], q)
    b0 = Poly(b.coeffs[:nhalf], q)
    b1 = Poly(b.coeffs[nhalf:], q)

    # Y = x^{n//2}
    # a = a_1 + Y a_0
    # b = b_1 + Y b_1
    # a*b = (a_1 + Y a_0)*(b_0 + Y b_1)
    #     = Y^2 a_1 b_1 + Y (a_0 b_1 + a_1 b_0) + a_0 b_0
    #     = Y^2 a_1 b_1 + Y ((a_0 + a_1)(b_0 + b_1) - a_1 b_1 - a_0 b_0) + a_0 b_0

    a0b0 = polymul_karatsuba_recursive(a0, b0, threshold)
    a1b1 = polymul_karatsuba_recursive(a1, b1, threshold)

    # compute t = ((a_0 + a_1)(b_0 + b_1) - a_1 b_1 - a_0 b_0)
    a0a1 = a0+a1
    b0b1 = b0+b1
    t = polymul_karatsuba_recursive(a0a1, b0b1, threshold)
    t -= a0b0
    t -= a1b1

    # compute a*b = Y^2 a_1 b_1 + Y t + a_0 b_0
    c = Poly.zero(2*n-1, q)
    for i in range(a0b0.n):
        c.coeffs[i] = a0b0.coeffs[i]

    for i in range(t.n):
        c.coeffs[i+nhalf] += t.coeffs[i]

    for i in range(a1b1.n):
        c.coeffs[i+2*nhalf] += a1b1.coeffs[i]

    # reduce coefficients mod q
    c.reduce()
    return c

def polymul_refined_karatsuba(a, b):
    """
        The core observation of refined Karatsuba is that some additions are performed twice.
        Consider the example with n =4.
        Let c = a0b0; d = (a0+a1)(b0+b1); e = a1b1;

        Normal Karatsuba computes:
          0    1    2    3    4    5    6
         ---- ---- ----      ---- ---- ----
        | c0 | c1 | c2 |    | e0 | e1 | e2 |
         ---- ---- ---- ---- ---- ----  ----
                + | d0 | d1 | d2 |
                   ---- ---- ----
                - | c0 | c1 | c2 |
                   ---- ---- ----
                - | e0 | e1 | e2 |
                   ---- ---- ----
        Here c2-e0 is computed twice (in column 2 and 4).
        For larger polynomials, this duplicate computation becomes significant (n//2-1)
        We can, thus, compute that part (say h) once and save some additions.
        Consequently, refined Karatsuba looks like
          0    1    2    3    4    5    6
         ---- ---- ----      ---- ---- ----
        | c0 | c1 | h0 |    |-h0 | e1 | e2 |
         ---- ---- ---- ---- ---- ----  ----
                + | d0 | d1 | d2 |
                   ---- ---- ----
                - | c0 | c1 |
                   ---- ---- ----
                -      | e1 | e2 |
                        ---- ----
    """

    assert a.n == b.n
    assert a.q == b.q
    n = a.n
    q = a.q

    nhalf = n//2

    a0 = Poly(a.coeffs[:nhalf], q)
    a1 = Poly(a.coeffs[nhalf:], q)
    b0 = Poly(b.coeffs[:nhalf], q)
    b1 = Poly(b.coeffs[nhalf:], q)

    a0b0 = a0*b0
    a1b1 = a1*b1

    # compute t = (a_0 + a_1)(b_0 + b_1)
    a0a1 = a0+a1
    b0b1 = b0+b1
    t = a0a1*b0b1

    # compute the sum that is used twice
    h  = Poly.zero(nhalf-1, q)
    for i in range(nhalf-1):
        # subtract lower half of a1b1 from upper half of a0b0
        h.coeffs[i] = (a0b0.coeffs[i+nhalf] - a1b1.coeffs[i]) % q

    # compute a*b = Y^2 a_1 b_1 + Y t + a_0 b_0
    c = Poly.zero(2*n-1, q)

    # in a real implementation, the copying can be mostly eliminated by using pointer magic.
    # copy lower part of a0b0
    for i in range(nhalf):
        c.coeffs[i] = a0b0.coeffs[i]

    # copy common sum
    for i in range(h.n):
        c.coeffs[nhalf+i] = h.coeffs[i]
        # this of course only helps us if we can do not have to do the negation
        # explicitly
        c.coeffs[2*nhalf+i] = -h.coeffs[i]

    # copy upper part of a1b1
    for i in range(a1b1.n - h.n):
        c.coeffs[i+2*nhalf+h.n] += a1b1.coeffs[i+h.n]

    # add entire (a_0 + a_1)(b_0 + b_1)
    for i in range(t.n):
        c.coeffs[i+nhalf] += t.coeffs[i]

    # subtract remaining part of a0b0
    for i in range(a0b0.n - h.n):
        c.coeffs[i+nhalf] -= a0b0.coeffs[i]

    # subtract remaining part of a1b1
    for i in range(a1b1.n - h.n):
        c.coeffs[i+nhalf+h.n] -= a1b1.coeffs[i+h.n]

    # reduce coefficients mod q
    c.reduce()
    return c


def testcase_karatsuba(n, q, printPoly=True):
    print(f"Testing Karatsuba with n={n}, q={q}")
    a = Poly.random(n, q)
    b = Poly.random(n, q)
    if printPoly: print("a=", a)
    if printPoly: print("b=", b)

    # compute reference product using schoolbook for comparison
    c_ref = a*b
    if printPoly: print("a*b (ref)=", c_ref)

    c = polymul_karatsuba(a, b)
    if printPoly: print("a*b (Karatsuba)=", c)
    print(f"equal: {c == c_ref}")
    if c == c_ref:
       return 0
    else:
       return 1

def testcase_karatsuba_recursive(n, q, t, printPoly=True):
    print(f"Testing recursive Karatsuba with n={n}, q={q}, t={4}")
    a = Poly.random(n, q)
    b = Poly.random(n, q)
    if printPoly: print("a=", a)
    if printPoly: print("b=", b)

    # compute reference product using schoolbook for comparison
    c_ref = a*b
    if printPoly: print("a*b (ref)=", c_ref)

    c = polymul_karatsuba_recursive(a, b, t)
    if printPoly: print("a*b (recursive Karatsuba)=", c)
    print(f"equal: {c == c_ref}")
    if c == c_ref:
       return 0
    else:
       return 1

def testcase_refined_karatsuba(n, q, printPoly=True):
    print(f"Testing refined Karatsuba with n={n}, q={q}")
    a = Poly.random(n, q)
    b = Poly.random(n, q)
    if printPoly: print("a=", a)
    if printPoly: print("b=", b)

    # compute reference product using schoolbook for comparison
    c_ref = a*b
    if printPoly: print("a*b (ref)=", c_ref)

    c = polymul_refined_karatsuba(a, b)
    if printPoly: print("a*b (refinedKaratsuba)=", c)
    print(f"equal: {c == c_ref}")
    if c == c_ref:
       return 0
    else:
       return 1

if __name__ == "__main__":
    rc = 0
    # plain Karatsuba tests
    rc |= testcase_karatsuba(n=5, q=17)
    rc |= testcase_karatsuba(n=256, q=1<<13, printPoly=False)

    # recursive Karatsuba tests
    rc |= testcase_karatsuba_recursive(n=8, q=17, t=2)
    rc |= testcase_karatsuba_recursive(n=256, q=1<<13, t=4, printPoly=False)

    # refined Karatsuba tests
    rc |= testcase_refined_karatsuba(n=4, q=17)
    rc |= testcase_refined_karatsuba(n=256, q=1<<13, printPoly=False)

    if rc != 0:
        print("TEST FAILED.")
        sys.exit(1)
    print("ALL GOOD.")