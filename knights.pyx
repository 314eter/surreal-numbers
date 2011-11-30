from nim cimport nimproduct, nimpower

cpdef unsigned knights(unsigned size, unsigned start, unsigned stop, unsigned previous):
    cdef unsigned n, t, test
    cdef bint bad
    for n in range(start, stop):
        bad = False
        for t in range(1, size):
            test = nimproduct(n, n ^ t)
            if test < size and test != previous:
                bad = True
                break
        if bad:
            continue
        if nimpower(n, size + 1) == previous:
            return n
