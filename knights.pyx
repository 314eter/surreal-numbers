#cython: cdivision=True
from nim cimport nimproduct
from libc.stdlib cimport malloc, free

cpdef unsigned knights(unsigned size, unsigned previous):
    cdef unsigned u, v, w, t, m, i, n0
    for u in range(log(size)):
        v = 2**u
        for w in range(v):
            t = v + w
            m = nimproduct(size / 2, nimproduct(t, t))
            for i in range(size / (2*v)):
                for n0 in range(2*v*i, 2*v*i + v):
                    if m ^ prod(n0, t) == previous:
                        return t * size + n0

prods = dict()

cdef unsigned prod(unsigned n, unsigned t):
    cdef unsigned res = 0
    cdef unsigned short bc = bitcount(n), exp = 0, count = 0
    cdef unsigned short* exps = <unsigned short*>malloc(bc * sizeof(unsigned short))
    while n:
        if n & 1:
            exps[count] = exp
            count += 1
        exp += 1
        n >>= 1
    for exp in exps[:bc]:
        if (t, exp) not in prods:
            prods[t, exp] = nimproduct(2**exp, 2**exp ^ t)
        res ^= prods[t, exp]
    free(exps)
    return res

cdef unsigned short log(unsigned n):
    cdef unsigned short l = 0
    while n:
        n <<= 1
        l += 1
    return l - 1

cdef unsigned short bitcount(unsigned n):
    cdef unsigned short count = 0
    while n:
        n &= n - 1
        count += 1
    return count
