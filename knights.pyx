#cython: cdivision=True
from nim cimport nimproduct
from libc.stdlib cimport malloc, free

cdef unsigned** prods

cdef unsigned short log(unsigned n):
    cdef unsigned short l = 0
    while n:
        n <<= 1
        l += 1
    return l - 1

cpdef unsigned knights(unsigned size, unsigned previous):
    cdef unsigned a, v, w, t, m, i, n0
    cdef unsigned short u, l = log(size)
    global prods
    prods = <unsigned**>malloc(size * sizeof(unsigned*))
    for a in range(size):
        prods[a] = <unsigned*>malloc(2 * l * sizeof(unsigned))
        for u in range(2 * l):
            prods[a][u] = 0
    for u in range(l):
        v = 2**u
        for w in range(v):
            t = v + w
            m = nimproduct(size / 2, nimproduct(t, t))
            for i in range(size / (2*v)):
                for n0 in range(2*v*i, 2*v*i + v):
                    if m ^ prod(n0, t) == previous:
                        for a in range(size):
                            free(prods[a])
                        free(prods)
                        return t * size + n0

cdef unsigned prod(unsigned n, unsigned t):
    cdef unsigned product, res = 0
    cdef unsigned short exp = 0
    global prods
    while n:
        if n & 1:
            if not prods[t][exp]:
                prods[t][exp] = nimproduct(2**exp, 2**exp ^ t)
            res ^= prods[t][exp]
        exp += 1
        n >>= 1
    return res
