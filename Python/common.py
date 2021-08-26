import math
from poly import Poly
from numpy import base_repr

# TODO: document each function
def bitrevidx(a, nbits):
    fmt = f"{{0:0{nbits}b}}"
    return list(map(lambda x: int(fmt.format(x)[::-1],2), a))

def bitreverse(a):
    b = [0]*len(a)
    logn =  int(math.log(len(a), 2))
    assert 2**logn == len(a)

    brv = bitrevidx(list(range(len(a))), logn)

    for i in range(len(a)):
        b[brv[i]] = a[i]
    return b

def revIdxBaseN(idx, ndigits, base):
    # represent as number in base N
    rep = base_repr(idx, base)
    # pad with zeros
    rep = "0"*(ndigits-len(rep))+ rep
    # reverse
    rep = rep[::-1]
    return int(rep, base)


def reverseBaseN(a, base):
    b = [0]*len(a)
    logn =  int(math.log(len(a), base))
    assert base**logn == len(a)

    for i in range(len(a)):
        b[revIdxBaseN(i, logn, base)] = a[i]
    return b

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

def ntt_naive_cyclic(p, root):
    q = p.q
    pntt = Poly.zero(p.n, q)
    for i in range(p.n):
        for j in range(p.n):
            pntt.coeffs[i] = (pntt.coeffs[i] +  p.coeffs[j]*pow(root, i*j, q)) % q
    return pntt

# TODO: describe what it is computing
def invntt_naive_cyclic(pntt, root):
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
    q = p.q
    pntt = Poly.zero(p.n, q)
    for i in range(p.n):
        for j in range(p.n):
            pntt.coeffs[i] = (pntt.coeffs[i] +  p.coeffs[j]*pow(root, j, q)*pow(root, 2*i*j, q)) % q
    return pntt

def invntt_naive_negacyclic(pntt, root):
    q = pntt.q
    p = Poly.zero(pntt.n, q)
    for i in range(pntt.n):
        for j in range(pntt.n):
            p.coeffs[i] = (p.coeffs[i] + pntt.coeffs[j]*pow(root, -2*i*j, q)) % q
    ninv = pow(p.n, -1, q)
    for i in range(p.n):
        p.coeffs[i] = (p.coeffs[i] * ninv * pow(root, -i, q)) % q
    return p


