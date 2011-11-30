cdef class Numbers(dict):
    """The Field of surreal numbers"""

    def __repr__(self):
        return 'Field of surreal numbers'

    def __call__(self, unsigned n):
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

    cpdef NumberIterator count(self, unsigned start=0, unsigned stop=0):
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

    def __add__(unsigned self, unsigned other):
        return N(self ^ other)

    def __sub__(unsigned self, unsigned other):
        return N(self ^ other)

    def __mul__(unsigned self, unsigned other):
        return N(nimproduct(self, other))

    def __div__(unsigned self, unsigned other):
        return N(nimproduct(self, niminvert(other)))

    def __pow__(self, int n, m):
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

    def __iter__(unsigned self):
        return NumberIterator(stop=self)

    def __len__(unsigned self):
        return self

    def __contains__(unsigned self, unsigned item):
        return item < self


cdef class NumberIterator:
    """Iterator over surreal numbers"""

    def __cinit__(self, unsigned start=0, unsigned stop=0):
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


cpdef inline unsigned nimsum(unsigned a, unsigned b):
    """Return the nim-sum of a and b."""
    return a ^ b

cdef object exps2(unsigned n):
    """Rewrite n as the sum of products of fermatpowers.
    
    A list of lists of exponents is returned such that:
    n = sum(2**(sum(2**n for n in exps)) for exps in exps2(n))
    
    """
    cdef unsigned short exp_n, exp_exp_n
    exponents = []
    exp_n = n.bit_length()
    for bit_n in bin(n)[2:]:
        exp_n -= 1
        if bit_n == '1':
            exponent = []
            exp_exp_n = exp_n.bit_length()
            for bit_exp_n in bin(exp_n)[2:]:
                exp_exp_n -= 1
                if bit_exp_n == '1':
                    exponent.append(exp_exp_n)
            exponents.append(exponent)
    return exponents

cpdef unsigned nimproduct(unsigned a, unsigned b):
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
    cdef unsigned result = 0
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

fermattable = {}

cpdef unsigned fermatproduct(object exps):
    """Return nim-product of one term (product of fermatpowers).

    Arguments:
    exps -- The exponents of the fermatpowers

    Algorithm:
    - Find least double exponent i in exps
    - If no doubles found, return 2**sum(2**n for n in exps)
    - Split exps in two list (a and b) without i, such that
      fermatproduct(exps) = fermatproduct(a) + fermatproduct(b)

    """
    cdef unsigned result = 0
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

cpdef unsigned nimpower(unsigned a, int n):
    """Return the power a**n based on nim-multiplication.

    Arguments:
    a -- base
    n -- exponent

    Algorithm:
    Exponentiation by squaring

    """
    cdef unsigned result = 1
    for b in bin(n)[2:]:
        result = nimproduct(result, result)
        if b == '1':
            result = nimproduct(result, a)
    return result

cpdef unsigned niminvert(unsigned a):
    """Return the nim-multiplicative inverse of a.

    Algorithm:
    As described in Ex. 5 of "Nim Multiplication" by H. W. Lenstra

    """
    cdef unsigned b = 2, x
    if a == 1:
        return 1
    while b**2 <= a:
        b = b**2
    x = a / b
    return nimproduct((a ^ x), niminvert(nimproduct(a, a ^ x)))
