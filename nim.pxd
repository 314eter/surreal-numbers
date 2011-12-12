cdef class NumberIterator:
    cdef unsigned long current
    cdef unsigned long stop

cdef class Numbers(dict):
    cpdef NumberIterator count(self, unsigned long start=*, unsigned long stop=*)

cdef class Number(int):
    cdef Number _inv
    cdef unsigned long _order

    cpdef unsigned long order(self)

cpdef inline unsigned long nimsum(unsigned long a, unsigned long b)

cpdef unsigned long nimproduct(unsigned long a, unsigned long b)

cpdef unsigned long nimpower(unsigned long a, long n)

cpdef unsigned long niminvert(unsigned long a)

cpdef unsigned long nimorder(unsigned long a)
