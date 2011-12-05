from nim cimport nimproduct

cpdef unsigned knights(unsigned size, unsigned previous,
                       unsigned start=1, unsigned stop=0):
    cdef unsigned t, m, n0
    if stop == 0:
        stop = size + 1
    checked = set()
    for t in range(start, stop):
        m = nimproduct(size / 2, nimproduct(t, t))
        checked.clear()
        for n0 in range(size):
            if n0 in checked:
                checked.remove(n0)
                continue
            checked.add(n0 ^ t)
            if m ^ nimproduct(n0, n0 ^ t) == previous:
                return t * size + n0

