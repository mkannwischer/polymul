""" Module containing common code for representing polynomials

"""

from enum import Enum
from random import randint

class Poly:
    """
    Representing polynomials in an arbitrary ring with coefficients modulo q


    Attributes
    ----------
    n : int
        number of coefficients (degree+1)
    coeffs : int[]
        the cofficients of the polynomial, with coeffs[i] corresponding to x^i
    """

    def __init__(self, coeffs, q=None):
        self.coeffs = coeffs
        self.n = len(coeffs)
        self.q = q

    def __str__(self):
        a = self.coeffs
        tmp = [f"{a[i]}x^{i}" for i in range(len(a)) if a[i] != 0][::-1]
        return " + ".join(tmp)


    def __mul__(self, other):
        """
        multiplies two polynomials and returns a double-sized product.
        Does not apply any modular reductions

        """
        if isinstance(other, Poly):
            assert self.q == other.q

            c = Poly.zero(self.n+other.n-1, self.q)
            for i in range(self.n):
                for j in range(other.n):
                    c.coeffs[i+j] += self.coeffs[i]*other.coeffs[j]

            # reduce coefficients modulo q
            c.reduce()
            return c
        elif isinstance(other, int):
            c = Poly.zero(self.n, self.q)
            for i in range(self.n):
                c.coeffs[i] = self.coeffs[i]*other
            c.reduce()
            return c
        else:
            raise NotImplementedError()

    def __rmul__(self, factor):
        return self*factor

    def __eq__(self, other):
        return self.coeffs == other.coeffs


    def __add__(self, other):
        assert self.q == other.q
        q = self.q
        c = Poly.zero(max(self.n, other.n), q)
        for i in range(self.n):
            c.coeffs[i] = self.coeffs[i]
        for i in range(other.n):
            c.coeffs[i] = c.coeffs[i]+other.coeffs[i]

        c.reduce()
        return c

    def __sub__(self, other):
        assert self.q == other.q
        q = self.q
        c = Poly.zero(max(self.n, other.n), q)
        for i in range(self.n):
            c.coeffs[i] = self.coeffs[i]
        for i in range(other.n):
            c.coeffs[i] = c.coeffs[i]-other.coeffs[i]

        c.reduce()
        return c

    def pointwise(self, other):
        assert self.q == other.q
        assert self.n == other.n
        c = Poly.zero(self.n, self.q)
        for i in range(self.n):
            c.coeffs[i] = self.coeffs[i] * other.coeffs[i]
        c.reduce()
        return c

    def __rshift__(self, shift):
        self.reduce()
        c = Poly.zero(self.n, self.q)
        for i in range(self.n):
            c.coeffs[i] = self.coeffs[i]>>shift
        return c

    def reduce(self):
        if self.q:
            for i in range(len(self.coeffs)):
                self.coeffs[i] %= self.q


    def copy(self):
        return Poly(self.coeffs.copy(), self.q)

    @staticmethod
    def random(n, q):
        """
        returns a uniformly random polynomial with n coefficients in [0,q)

        Parameters
        ----------
        n : int
            number of coefficients
        q : int
            modulus
        """

        return Poly([randint(0, q-1) for _ in range(n)], q)

    @staticmethod
    def zero(n, q):
        """
        returns a polynomial with n all-zero coefficients

        Parameters
        ----------
        n : int
            number of coefficients
        """
        return Poly([0 for _ in range(n)], q)

