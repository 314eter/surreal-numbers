#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: overflowcheck=False
from nim cimport nimproduct
from libc.stdlib cimport malloc, free
from libc.string cimport memset

cdef unsigned long* prods

cpdef unsigned long knight(unsigned short n, unsigned long previous, unsigned long start=0):
    cdef unsigned long size, v, w, t, m, i, n0
    cdef unsigned short exp, u, l = 2**n
    global prods
    prods = <unsigned long*>malloc(2 * l * sizeof(unsigned long))
    size = 2**l
    for u in range(l / 4, l):
        v = 2**u
        for w in range(v):
            t = v + w
            print 't:', t
            memset(prods, 0, 2 * l * sizeof(unsigned long))
            m = nimproduct(size >> 1, nimproduct(t, t))
            for i in range(size / (2*v)):
                for n0 in range(2*v*i, 2*v*i + v):
                    if m ^ prod(n0, t) == previous:
                        free(prods)
                        return t * size + n0

cdef unsigned long prod(unsigned long n, unsigned long t):
    cdef unsigned long res = 0
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
