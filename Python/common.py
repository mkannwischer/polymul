import math
from poly import Poly
from numpy import base_repr



def isPrime(n):
    """
        checks if a number is prime

        Parameters
        ----------
        n : int
    """

    for i in range(2, n):
        if n%i == 0:
            return False
    return True

def primitiveRootOfUnity(n, q):
        """
            returns a primitive n-th root of unity mod q if it exists, i.e.,
            returns w, s.t. w^n = 1 mod q and w^k != 1 mod q for all k < n.
            Raises an error if it does not exist.


        """

    for i in range(2,q):
        if pow(i, n, q) == 1:
            # i is an n-th root of unity, but it may not be primitive
            # check by making sure i^j != 1 mod q for 1<=j<n
            isPrimitive = True
            for j in range(1,n):
                if pow(i, j, q) == 1:
                    isPrimitive = False
                    break
            if isPrimitive:
                return i
    raise Exception(f"{n}th root of unity mod {q} does not exist")

def bitreverse(a):
        """
            Transforms an array a into bitreversed order.
            a needs to have length 2^k for some k
            For example
            [0,1,2,3,4,5,6,7] will be turned into [0,4,2,6,1,5,3,7]

        """

    b = [0]*len(a)
    logn =  int(math.log(len(a), 2))
    assert 2**logn == len(a)

    def bitrevidx(a, nbits):
        fmt = f"{{0:0{nbits}b}}"
        return list(map(lambda x: int(fmt.format(x)[::-1],2), a))

    brv = bitrevidx(list(range(len(a))), logn)

    for i in range(len(a)):
        b[brv[i]] = a[i]
    return b


def reverseBaseN(a, base):
        """
            Transforms an array into reversed order for any base.
            For base=2, this is the same as bitreverse.
            a needs to have length base^k for some k.
            For example
            [0,1,2,3,4,5,6,7,8] for base=3 will be turned into [0,3,6,1,4,7,2,5,8]

        """

    b = [0]*len(a)
    logn =  int(math.log(len(a), base))
    assert base**logn == len(a)

    def revIdxBaseN(idx, ndigits, base):
    # represent as number in base N
        rep = base_repr(idx, base)
        # pad with zeros
        rep = "0"*(ndigits-len(rep))+ rep
        # reverse
        rep = rep[::-1]
        return int(rep, base)

    for i in range(len(a)):
        b[revIdxBaseN(i, logn, base)] = a[i]
    return b



def ntt_naive_cyclic(p, root):
        """
            Naive implementation of a cyclic NTT. Needs an n-th root of unity.
            Computes
            antt_i = sum_j=0^n (a_j  root^(ij)) for 0 <= i < n

        """
    
    q = p.q
    pntt = Poly.zero(p.n, q)
    for i in range(p.n):
        for j in range(p.n):
            pntt.coeffs[i] = (pntt.coeffs[i] +  p.coeffs[j]*pow(root, i*j, q)) % q
    return pntt

def invntt_naive_cyclic(pntt, root):
        """
            Naive implementation of an inverse cyclic NTT. Needs an n-th root of unity.
            Computes
            a_i = 1/n * sum_j=0^n (a_j  root^(-ij)) for 0 <= i < n

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

def ntt_naive_negacyclic(p, root):
        """
            Naive implementation of a negacyclic NTT. Needs a 2n-th root of unity.
            Computes
            antt_i = sum_j=0^n (a_j  root^(2ij + j)) for 0 <= i < n

        """
    q = p.q
    pntt = Poly.zero(p.n, q)
    for i in range(p.n):
        for j in range(p.n):
            pntt.coeffs[i] = (pntt.coeffs[i] +  p.coeffs[j]*pow(root, j, q)*pow(root, 2*i*j, q)) % q
    return pntt

def invntt_naive_negacyclic(pntt, root):
        """
            Naive implementation of an inverse negacyclic NTT. Needs a 2n-th root of unity.
            Computes
            a_i = 1/n * root^(-i) * sum_j=0^n (a_j  root^(-2ij)) for 0 <= i < n

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


