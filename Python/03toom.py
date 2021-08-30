#!/usr/bin/env python3
"""Section 2.2.3: Toom--Cook multiplication.

Illustrates Toom-3, and Toom-4 multiplication for any modulus.
"""

from poly import Poly
import sys

def polymul_toom3(a, b):
    """Compute polynomial product using three-way Toom--Cook multiplication.

    Assumes n is the same for both a and b; also 3 must divide n.
    Sets y=x^(n/3), s.t.
    a = a0 + y a1 + y^2 a2
    b = b0 + y b1 + y^2 b2

    Then evaluates at y = {0, infty, 1, -1, -2}, performs 5 multiplications
    of degree n/3 polynomials, and interpolates c using the formulas presented
    in the thesis.

    During the interpolation, we need to divide by 2 and 3
    In case q is co-prime to 2 and 3, we can simply multiply by the inverse

    As q is often a power of two, the inverse of two does often not exist.
    In that case there is a workaround: We do the small multiplications
    modulo 2q instead of q, such that we can divide by two during the
    interpolation without changing the result.


    Parameters
    ----------
    a : Poly
        first multiplicand with n coefficients.
    b : Poly
        second multiplicand with n coefficients.
    Returns
    ----------
    Poly
        product a*b with 2n-1 coefficients
    """
    assert a.n == b.n
    assert a.q == b.q
    q = a.q
    n = a.n
    # for simplicity, we assume that inputs nicely splits into three parts of same size
    assert n % 3 == 0
    # we also assume that q is co-prime to 3, s.t., the inverse of 3 exists mod q
    assert q % 3 != 0

    qq = q
    try:
        inv2 = pow(2, -1, q)
    except ValueError:
        inv2 = None
        qq = 2*q

    inv3 = pow(3, -1, qq)

    # split a and b in 3 parts each, i.e.,
    # a = a0 + y a1 + y^2 a2
    # b = b0 + y b1 + y^2 b2
    nsplit = n//3
    a0 = Poly(a.coeffs[:nsplit], qq)
    a1 = Poly(a.coeffs[nsplit:2*nsplit], qq)
    a2 = Poly(a.coeffs[2*nsplit:], qq)
    b0 = Poly(b.coeffs[:nsplit], qq)
    b1 = Poly(b.coeffs[nsplit:2*nsplit], qq)
    b2 = Poly(b.coeffs[2*nsplit:], qq)

    # Evaluate at y={0, \infty, 1, -1, -2}
    a_0   = a0
    a_inf = a2
    a_1   = a0 + a1 + a2
    a_m1  = a0 - a1 + a2
    a_m2  = a0 - 2*a1 + 4*a2

    b_0   = b0
    b_inf = b2
    b_1   = b0 + b1 + b2
    b_m1  = b0 - b1 + b2
    b_m2  = b0 - 2*b1 + 4*b2

    # do smaller multiplications (in this case using schoolbook; could also be
    # other Toom or Karatsuba)
    c_0   = a_0 * b_0
    c_inf = a_inf * b_inf
    c_1   = a_1 * b_1
    c_m1  = a_m1 * b_m1
    c_m2  = a_m2 * b_m2

    # interpolation of c = c_0 + v1  y^1 + v2 y^2 + v3 y^3 + c_inf y^4
    if inv2 != None:
        # if 2 is invertible, we multiply by the inverse
        v3 = (c_m2 - c_1) * inv3            # -v1 + v2 - 3v3 + 5c_inf
        v1 = (c_1 - c_m1) * inv2            #  v1      +  v3
        v2 = c_m1 - c_0                     # -v1 + v2 -  v3 +  c_inf
        v3 = (v2 - v3) * inv2 + 2 * c_inf   #             v3
        v2 = v2 + v1 - c_inf                #       v2
        v1 = v1 - v3                        #  v1
    else:
        # otherwise, we divide by two (in which case the small multiplications need
        # to be correct mod 2*q)
        v3 = (c_m2 - c_1) * inv3            # -v1 + v2 - 3v3 + 5c_inf
        v1 = (c_1 - c_m1) >> 1              #  v1      +  v3
        v2 = c_m1 - c_0                     # -v1 + v2 -  v3 +  c_inf
        v3 = ((v2 - v3)>>1) + 2 * c_inf     #             v3
        v2 = v2 + v1 - c_inf                #       v2
        v1 = v1 - v3                        #  v1


    # recombine 5 (overlapping) parts
    # switch back to a small modulus q now
    c = Poly.zero(n*2-1, q)
    for i in range(c_inf.n):
        c.coeffs[0*nsplit+i] += c_0.coeffs[i]
        c.coeffs[1*nsplit+i] += v1.coeffs[i]
        c.coeffs[2*nsplit+i] += v2.coeffs[i]
        c.coeffs[3*nsplit+i] += v3.coeffs[i]
        c.coeffs[4*nsplit+i] += c_inf.coeffs[i]
    c.reduce()
    return c

def polymul_toom4(a, b):
    """Compute polynomial product using four-way Toom--Cook multiplication.

    Assumes n is the same for both a and b; also 4 must divide n.
    Sets y=x^(n/4), s.t.
    a = a0 + y a1 + y^2 a2 + y^3 a3
    b = b0 + y b1 + y^2 b2 + y^3 a3

    Then evaluates at y = {0, infty, 1, -1, 2, -2, 3}, performs 7 multiplications
    of degree n/4 polynomials, and interpolates c using the formulas presented
    in the thesis.

    During the interpolation, we need to divide by 8 and 3
    In case q is co-prime to 8 and 3, we can simply multiply by the inverse.
    As q is often a power of two, the inverse of 8 does often not exist.
    In that case there is a workaround: We do the small multiplications
    modulo 8q instead of q, such that we can divide by two during the
    interpolation without changing the result.

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
    q = a.q
    n = a.n
    # for simplicity, we assume that inputs nicely splits into four parts of same size
    assert n % 4 == 0
    assert q % 3 != 0

    qq = q
    try:
        inv8 = pow(8, -1, q)
        inv2 = pow(2, -1, q)
        inv4 = pow(4, -1, q)
    except ValueError:
        inv8 = inv4 = inv2 = None
        qq = 8*q

    inv3 = pow(3, -1, qq)
    inv5 = pow(5, -1, qq)

    # split a and b in 4 parts each, i.e.,
    # a = a0 + y a1 + y^2 a2 + y^3 a3
    # b = b0 + y b1 + y^2 b2 + y^3 b3
    nsplit = n//4
    a0 = Poly(a.coeffs[:nsplit], qq)
    a1 = Poly(a.coeffs[nsplit:2*nsplit], qq)
    a2 = Poly(a.coeffs[2*nsplit:3*nsplit], qq)
    a3 = Poly(a.coeffs[3*nsplit:], qq)
    b0 = Poly(b.coeffs[:nsplit], qq)
    b1 = Poly(b.coeffs[nsplit:2*nsplit], qq)
    b2 = Poly(b.coeffs[2*nsplit:3*nsplit], qq)
    b3 = Poly(b.coeffs[3*nsplit:], qq)

    # Evaluate at y={0, \infty, 1, -1, 2, -2, 3}
    a_0  = a0
    a_inf= a3
    a_1  = a0 +   a1 +   a2 +    a3
    a_m1 = a0 -   a1 +   a2 -    a3
    a_2  = a0 + 2*a1 + 4*a2 +  8*a3
    a_m2 = a0 - 2*a1 + 4*a2 -  8*a3
    a_3  = a0 + 3*a1 + 9*a2 + 27*a3

    b_0  = b0
    b_inf= b3
    b_1  = b0 +   b1 +   b2 +    b3
    b_m1 = b0 -   b1 +   b2 -    b3
    b_2  = b0 + 2*b1 + 4*b2 +  8*b3
    b_m2 = b0 - 2*b1 + 4*b2 -  8*b3
    b_3  = b0 + 3*b1 + 9*b2 + 27*b3

    # do smaller multiplications (in this case using schoolbook; could also be
    # other Toom or Karatsuba)
    c_0   = a_0 * b_0
    c_inf = a_inf * b_inf
    c_1   = a_1 * b_1
    c_m1  = a_m1 * b_m1
    c_2   = a_2 * b_2
    c_m2  = a_m2 * b_m2
    c_3   = a_3 * b_3

    # Interpolation of c = c_0 + v1  y^1 + v2  y^2 + v3  y^3 + v4  y^4 + v5 y^5
    #                     + c_inf y^6
    if inv8 != None:
        # if 8,4,2 is invertible, we multiply by the inverses
        t0 = (c_1 + c_m1) * inv2 - c_0 - c_inf            # v2 +  v4
        t1 = (c_2 + c_m2 - 2*c_0 - 128*c_inf) * inv8      # v2 + 4v4
        v4 = (t1 - t0)*inv3                               # v4
        v2 = (t0 - v4)                                    # v2
        t0 = (c_1 - c_m1) * inv2                          # v1 + v3 + v5
        t1 = ((c_2 - c_m2)*inv4 - t0)* inv3               # v3 + 5v5
        t2 = (c_3 - c_0 - 9*v2 - 81*v4 - 729*c_inf)*inv3  # v1 + 9v3 + 81v5
        t2 = (t2 - t0)*inv8 - t1                          # 5v5
        v5 = t2*inv5                                      # v5
        v3 = t1 - t2                                      # v3
        v1 = t0 - v3 - v5                                 # v1
    else:
        # otherwise, we divide (in which case the small multiplications need
        # to be correct mod 8*q; and the interpolation needs to be done mod 16*q)
        c_0.q = c_1.q = c_m1.q = c_2.q = c_m2.q = c_3.q = c_inf.q = 2*qq

        t0 = ((c_1 + c_m1)>>1) - c_0 - c_inf              # v2 +  v4
        t1 = (c_2 + c_m2 - 2*c_0 - 128*c_inf)>>3          # v2 + 4v4
        v4 = (t1 - t0)*inv3                               # v4
        v2 = (t0 - v4)                                    # v2
        t0 = (c_1 - c_m1)>>1                              # v1 + v3 + v5
        t1 = (((c_2 - c_m2)>>2) - t0)* inv3               # v3 + 5v5
        t2 = (c_3 - c_0 - 9*v2 - 81*v4 - 729*c_inf)*inv3  # v1 + 9v3 + 81v5
        t2 = ((t2 - t0)>>3) - t1                          # 5v5
        v5 = t2*inv5                                      # v5
        v3 = t1 - t2                                      # v3
        v1 = t0 - v3 - v5                                 # v1

    # recombine 7 (overlapping) parts
    # switch back to a small modulus q
    c = Poly.zero(n*2-1, q)
    for i in range(c_inf.n):
        c.coeffs[0*nsplit+i] += c_0.coeffs[i]
        c.coeffs[1*nsplit+i] += v1.coeffs[i]
        c.coeffs[2*nsplit+i] += v2.coeffs[i]
        c.coeffs[3*nsplit+i] += v3.coeffs[i]
        c.coeffs[4*nsplit+i] += v4.coeffs[i]
        c.coeffs[5*nsplit+i] += v5.coeffs[i]
        c.coeffs[6*nsplit+i] += c_inf.coeffs[i]
    c.reduce()
    return c

def testcase_toom3(n, q, printPoly=True):
    """Random test of three-way Toom--Cook multiplication (`polymul_toom3`).

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
    print(f"Testing Toom3 multiplication with q={q}, n={n}")
    a = Poly.random(n, q)
    b = Poly.random(n, q)
    if printPoly: print("a=", a)
    if printPoly: print("b=", b)

    c_ref = a*b
    if printPoly: print("a*b (ref)=", c_ref)

    c = polymul_toom3(a, b)
    if printPoly: print("a*b (toom3)=", c)
    print(f"equal: {c == c_ref}")
    if c == c_ref:
        return 0
    else:
        return 1

def testcase_toom4(n, q, printPoly=True):
    """Random test of four-way Toom--Cook multiplication (`polymul_toom4`).

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
    print(f"Testing Toom4 multiplication with q={q}, n={n}")
    a = Poly.random(n, q)
    b = Poly.random(n, q)
    if printPoly: print("a=", a)
    if printPoly: print("b=", b)

    c_ref = a*b
    if printPoly: print("a*b (ref)=", c_ref)

    c = polymul_toom4(a, b)
    if printPoly: print("a*b (toom4)=", c)
    print(f"equal: {c == c_ref}")
    if c == c_ref:
        return 0
    else:
        return 1


if __name__ == "__main__":
    rc = 0

    # toom 3
    rc |= testcase_toom3(n=6, q=17)                        # inverses of 2 and 3 exist
    rc |= testcase_toom3(n=6, q=8)                         # inverse of 2 doesn't exist
    rc |= testcase_toom3(n=258, q=1<<13, printPoly=False)  # inverse of 2 doesn't exist

    # toom 4`
    rc |= testcase_toom4(n=8, q=17)                        # inverses of 2,3,4,5,8 exist
    rc |= testcase_toom4(n=8, q=1<<13)                     # inverses of 2,4,8 don't exist
    rc |= testcase_toom4(n=256, q=17, printPoly=False)     # inverses of 2,3,4,5,8 exist
    rc |= testcase_toom4(n=256, q=1<<13, printPoly=False)  # inverses of 2,3,4,5,8 exist

    if rc != 0:
        print("TEST FAILED.")
        sys.exit(1)
    print("ALL GOOD.")