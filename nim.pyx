#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: overflowcheck=False
cdef class NumberIterator:
    """Iterator over surreal numbers"""

    def __cinit__(self, unsigned long start=0, unsigned long stop=0):
        """Iterate over the interval [start, stop[."""
        self.current = start
        self.stop = stop

    def __iter__(self):
        return self

    def __next__(self):
        if self.stop and self.current == self.stop:
            raise StopIteration
        self.current += 1
        return N(self.current - 1)


cdef class Numbers(dict):
    """The Field of surreal numbers"""

    def __repr__(self):
        return 'Field of surreal numbers'

    def __call__(self, unsigned long n):
        """Return the number n

        Arguments:
        n -- an integer to convert to a surreal number

        """
        if n not in self:
            self[n] = Number(n)
        return self[n]

    def __iter__(self):
        """Return an iterator over all numbers."""
        return NumberIterator()

    cpdef NumberIterator count(self, unsigned long start=0, unsigned long stop=0):
        """Return an iterator over a range of Numbers.

        Arguments:
        start -- the first number (default 0)
        stop -- the first number not included (default 0)

        """
        return NumberIterator(start, stop)

N = Numbers()


cdef class Number(int):
    """A surreal number

    This class inherits from int, and overrides +-*/ with nim-operations.
    The preferred way to initialize is N(a).

    Operations:
    a + b -- nim-addition
    a - b -- nim-addition
    a * b -- nim-multiplication
    a / b -- nim-multiplication of a and ~b
    a**n  -- exponentiation based on nim-multiplication
    ~a    -- invers (a * ~a = 1)

    Iteration over a iterates over all numbers less than a

    """

    def __repr__(self):
        return 'N(' + int.__repr__(self) + ')'

    def __add__(unsigned long self, unsigned long other):
        return N(self ^ other)

    def __sub__(unsigned long self, unsigned long other):
        return N(self ^ other)

    def __mul__(unsigned long self, unsigned long other):
        return N(nimproduct(self, other))

    def __div__(unsigned long self, unsigned long other):
        return N(nimproduct(self, niminvert(other)))

    def __pow__(self, long n, m):
        if n < 0:
            self = ~self
            n = -n
        return N(nimpower(self, n))

    def __neg__(self):
        return self

    def __pos__(self):
        return self

    def __invert__(self):
        if self == 0:
            raise ZeroDivisionError
        if not self._inv:
            self._inv = N(niminvert(self))
        return self._inv

    def __iter__(unsigned long self):
        return NumberIterator(stop=self)

    def __len__(unsigned long self):
        return self

    def __contains__(unsigned long self, unsigned long item):
        return item < self

    cpdef unsigned long order(self):
        if self == 0:
            raise ZeroDivisionError
        if not self._order:
            self._order = nimorder(self)
        return self._order


cpdef inline unsigned long nimsum(unsigned long a, unsigned long b):
    """Return the nim-sum of a and b."""
    return a ^ b

cdef dict expstable = dict()

cdef object exps2(unsigned long n):
    """Rewrite n as the sum of products of fermatpowers.

    A list of lists of exponents is returned such that:
    n = sum(2**(sum(2**n for n in exps)) for exps in exps2(n))

    """
    cdef unsigned short i = 0, j
    cdef unsigned long m = 1
    if n in expstable:
        return expstable[n]
    exponents = []
    while n >= m:
        if n & m:
            exponent = []
            j = 0
            while i >= (1 << j):
                if i & (1 << j):
                    exponent.append(j)
                j += 1
            exponents.append(exponent)
        i += 1
        m = m << 1
    expstable[n] = exponents
    return exponents

cpdef unsigned long nimproduct(unsigned long a, unsigned long b):
    """Return the nim-product of a and b.

    Algorithm:
    - Rewrite a and b with exps2
    - Compute product with distributivuty of multiplication over addition.
      The resulting list of list of exponents f has the same property as exps2:
      a * b = sum(2**(sum(2**n for n in exps)) for exps in f)
    - Remove double terms (because a + a = 0)
    - Compute nim-product of each term (product of fermatpowers)
    - Compute nim-sum of terms

    """
    cdef unsigned long result = 0
    f = set()
    for l in exps2(a):
        for r in exps2(b):
            exps = tuple(sorted(l + r, reverse=True))
            if exps in f:
                f.remove(exps)
            else:
                f.add(exps)
    for exps in f:
        result ^= fermatproduct(exps)
    return result

cdef dict fermattable = dict()

cdef unsigned long fermatproduct(object exps):
    """Return nim-product of one term (product of fermatpowers).

    Arguments:
    exps -- The exponents of the fermatpowers

    Algorithm:
    - Find least double exponent i in exps
    - If no doubles found, return 2**sum(2**n for n in exps)
    - Split exps in two list (a and b) without i, such that
      fermatproduct(exps) = fermatproduct(a) + fermatproduct(b)

    """
    cdef unsigned long result = 0
    cdef int i = len(exps) - 1
    while i > 0 and exps[i] != exps[i - 1]:
        i -= 1
    if i <= 0:
        for i in exps:
            result += 2**i
        return 2**result
    if exps in fermattable:
        return fermattable[exps]
    result = fermatproduct(exps[:i-1] + exps[i:])
    result ^= fermatproduct(tuple(sorted(list(exps[:i-1] + exps[i+1:]) +
                            range(exps[i]-1, -1, -1), reverse=True)))
    fermattable[exps] = result
    return result

cpdef unsigned long nimpower(unsigned long a, long n):
    """Return the power a**n based on nim-multiplication.

    Arguments:
    a -- base
    n -- exponent

    Algorithm:
    Exponentiation by squaring

    """
    cdef unsigned long i = 0, result = 1
    for b in bin(n)[2:]:
        result = nimproduct(result, result)
        if b == '1':
            result = nimproduct(result, a)
    return result

cpdef unsigned long niminvert(unsigned long a):
    """Return the nim-multiplicative inverse of a.

    Algorithm:
    As described in Ex. 5 of "Nim Multiplication" by H. W. Lenstra

    """
    cdef unsigned long x, b = 2
    if a == 1:
        return 1
    while b**2 <= a:
        b = b**2
    x = a / b
    return nimproduct((a ^ x), niminvert(nimproduct(a, a ^ x)))

cpdef unsigned long nimorder(unsigned long a):
    cdef unsigned long exp = 0
    cdef unsigned long field
    if a == 1:
        return 1
    while 2**(2**exp) <= a:
        exp += 1
    field = 2**(2**exp) - 1
    for exp in range(3, field):
        if field % exp == 0 and nimpower(a, exp) == 1:
            return exp
    return field
