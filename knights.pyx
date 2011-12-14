#cython: cdivision=True
from nim cimport nimproduct
from libc.stdlib cimport calloc, free

cdef unsigned long* prods

cdef unsigned short log(unsigned long n):
    cdef unsigned short l = 0
    while n:
        n <<= 1
        l += 1
    return l - 1

cpdef unsigned long knights(unsigned long size, unsigned long previous):
    cdef unsigned long v, w, t, m, i, n0
    cdef unsigned short exp, u, l = log(size)
    global prods
    for u in range(l):
        v = 2**u
        for w in range(v):
            t = v + w
            prods = <unsigned long*>calloc(2 * l, sizeof(unsigned long))
            m = nimproduct(size / 2, nimproduct(t, t))
            for i in range(size / (2*v)):
                for n0 in range(2*v*i, 2*v*i + v):
                    if m ^ prod(n0, t) == previous:
                        free(prods)
                        return t * size + n0
            free(prods)

cdef unsigned long prod(unsigned long n, unsigned long t):
    cdef unsigned long product, res = 0
    cdef unsigned short exp = 0
    global prods
    while n:
        if n & 1:
            if not prods[exp]:
                prods[exp] = nimproduct(2**exp, 2**exp ^ t)
            res ^= prods[exp]
        exp += 1
        n >>= 1
    return res
