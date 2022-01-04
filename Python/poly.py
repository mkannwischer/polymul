"""Module containing common code for representing polynomials."""

from enum import Enum
from random import randint

class Poly:
    """Representing polynomials in an arbitrary ring with coefficients modulo q.

    Attributes
    ----------
    n : int
        number of coefficients (degree+1).
    coeffs : int[]
        the cofficients of the polynomial, with coeffs[i] corresponding to x^i.
    q : int
        modulus. Optional.
    """

    def __init__(self, coeffs, q=None):
        """Initialize Poly given a list of coefficients and optionally a modulus q."""
        self.coeffs = coeffs
        self.n = len(coeffs)
        self.q = q

    def __str__(self):
        """Convert a Poly to a string e.g., '1x^2+2x+3'."""
        a = self.coeffs
        tmp = [f"{a[i]}x^{i}" for i in range(len(a)) if a[i] != 0][::-1]
        return " + ".join(tmp)


    def __mul__(self, other):
        """Multiply two polynomials and returns a double-sized product or multiplies a polynomial by an integer.

        Does not apply any modular reductions mod q
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
        """multiply a polynomial by an integer."""
        return self*factor

    def __eq__(self, other):
        """Check two polynomials for equality by comparing their coefficient lists."""
        return self.coeffs == other.coeffs


    def __add__(self, other):
        """Add two polynomials and perform a reduction modulo q."""
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
        """Subtract two polynomials and perform a reduction modulo q."""
        assert self.q == other.q
        q = self.q
        c = Poly.zero(max(self.n, other.n), q)
        for i in range(self.n):
            c.coeffs[i] = self.coeffs[i]
        for i in range(other.n):
            c.coeffs[i] = c.coeffs[i]-other.coeffs[i]

        c.reduce()
        return c

    def __rshift__(self, shift):
        """Shift each coefficient of a polynomial to the right."""
        self.reduce()
        c = Poly.zero(self.n, self.q)
        for i in range(self.n):
            c.coeffs[i] = self.coeffs[i]>>shift
        return c


    def pointwise(self, other):
        """Perform a coefficient-wise multiplication followed by a reduction modulo q.

        Parameters
        ----------
        other: Poly
               The other polynomial.

        """
        assert self.q == other.q
        assert self.n == other.n
        c = Poly.zero(self.n, self.q)
        for i in range(self.n):
            c.coeffs[i] = self.coeffs[i] * other.coeffs[i]
        c.reduce()
        return c


    def reduce(self):
        """Reduce a polynomial mod q, i.e., [0,q)."""
        if self.q:
            for i in range(len(self.coeffs)):
                self.coeffs[i] %= self.q


    def copy(self):
        """Create a copy of a polynomial."""
        return Poly(self.coeffs.copy(), self.q)

    @staticmethod
    def random(n, q):
        """Sample a uniformly random polynomial with n coefficients in [0,q).

        Parameters
        ----------
        n : int
            number of coefficients.
        q : int
            modulus.
        """
        return Poly([randint(0, q-1) for _ in range(n)], q)

    @staticmethod
    def zero(n, q):
        """Create a polynomial with coefficients all n coefficients zero.

        Parameters
        ----------
        n : int
            number of coefficients.
        """
        return Poly([0 for _ in range(n)], q)

