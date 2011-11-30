cdef class Numbers(dict):
    cpdef NumberIterator count(self, unsigned start=*, unsigned stop=*)

cdef class Number(int):
    cdef Number _inv

cdef class NumberIterator:
    cdef unsigned current
    cdef unsigned stop

cpdef inline unsigned nimsum(unsigned a, unsigned b)

cpdef unsigned nimproduct(unsigned a, unsigned b)

cpdef unsigned nimpower(unsigned a, int n)

cpdef unsigned niminvert(unsigned a)
