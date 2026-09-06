from PyPR.FeedbackFunctions import FeedbackFunction
from PyPR.BooleanLogic import BooleanFunction, BooleanANF

from PyPR import FeedbackRegister
from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey

from functools import cached_property

class MPR(FeedbackFunction):
    """A Mersenne Product Register over the finite field GF(2^n).

    An MPR identifies its `size`-bit state with an element of the field
    GF(2^n) = GF(2)[x] / P(x), where P is a primitive polynomial of degree n.
    Each clock cycle multiplies the current state by a fixed update polynomial
    U modulo P. With the default U(x) = x (multiplication by the field
    generator) the orbit of any nonzero state has length 2^n - 1 (a Mersenne
    number when 2^n - 1 is prime, hence the name). Other choices of U produce
    a register that traces the same cycle but in a permuted order, altering
    the algebraic relationships between successive states without changing
    the period.

    Conceptually this is a "shift up" register: multiplication by x sends
    α^k to α^{k+1}, the opposite direction from Galois and Fibonacci LFSRs in
    this library, which effectively multiply by x^{-1} (shift down). In the
    primal/dual disambiguation used elsewhere in the library, the polynomial
    P passed to an MPR is the primal — the output sequence satisfies P
    directly with no reversal — whereas Galois and Fibonacci take their
    input as the dual. See `docs/conventions/Polynomial Conventions.md` for the
    precise distinction.

    :ivar size: The number of bits in the register (the field extension degree n).
    :vartype size: int
    :ivar primitive_polynomial: The coefficient list of P, of length n+1,
        with `primitive_polynomial[i]` the coefficient of x^i.
    :vartype primitive_polynomial: list[int]
    :ivar update_polynomial: The coefficient list of U, of length n, with
        `update_polynomial[i]` the coefficient of x^i.
    :vartype update_polynomial: list[int]
    """

    primitive_polynomial: list[int]
    update_polynomial: list[int]

    def __init__(self,
        size: int,
        primitive_poly: str | list[int],
        update_poly: list[int] | None = None
    ) -> None:
        """Construct an MPR of the given size with the specified primitive
        polynomial and update polynomial.

        The state evolves under multiplication by U (modulo P) in the field
        GF(2^n) = GF(2)[x] / P(x). Each output bit function is built
        symbolically by first carrying out the polynomial multiplication
        U * state in GF(2)[x] (which can push terms up to degree 2n - 2), then
        reducing modulo P using the relation x^n = P(x) - x^n to express x^n
        through x^{2n-2} as linear combinations of x^0 through x^{n-1}.

        :param size: The number of bits in the register, equal to the field
            extension degree n.
        :type size: int
        :param primitive_poly: The primitive polynomial P of degree n (the
            field modulus). May be given either as a coefficient list of
            length n+1 (with index i holding the coefficient of x^i, so the
            leading coefficient at index n is 1) or as a Koopman hex string.
            The Koopman format encodes the polynomial's lower n coefficients
            in hex; the leading 1 at degree n is implicit (e.g. "12" ->
            binary "10010" -> "100101" = 1 + x^3 + x^5).
        :type primitive_poly: str | list[int]
        :param update_poly: The update polynomial U of degree less than n,
            given as a coefficient list of length n. Defaults to U(x) = x
            (multiplication by the field generator), under which the MPR
            realizes the standard cyclic group action on nonzero field
            elements.
        :type update_poly: list[int] | None
        """
        self.size = size
        self._data_store = None

        # Format Update & Primitive polyomials: -----------------------------

        # convert U to a list:
        if not update_poly:
            update_poly = [0,1] + [0]*(size-2)
        elif type(update_poly) == int:
            update_poly = [int(x) for x in format(int(primitive_poly,16), f"0>{size}b")]

        # convert update to powers ([1,0,0,1,1] -> [0,3,4])
        self.update_polynomial = update_poly
        update_powers = [idx for (idx, t) in enumerate(update_poly) if t == 1]

        # P can be either polynomial list or koopman hex string:
        #   -koopman format note: binary interpretation is missing final 1
        #   -example: "12" -> "10010" -> "100101" = (1 + x^3 + x^5) -> [0,3,5]
        if isinstance(primitive_poly, str):
            primitive_poly = [int(x) for x in format(int(primitive_poly,16), f"0>{size}b")] + [1]
            # primitive_poly = list(BitVector(intVal = int(primitive_poly, 16), size = size))+[1]
            
        self.primitive_polynomial = primitive_poly
        primitive_powers = [(idx) for (idx, t) in enumerate(primitive_poly) if t == 1]

        # symbolically represent multiplication in GF(2^n): ----------------------------

        # the anf also includes n-1 "hypothetical bits" for higher powers
        # anf[:size] are real bits, anf [size:] are hypothetical bits
        functions = [[] for i in range(2*size)]
        
        # Multiply by update polynomial U
        for idx in range(size):
            for power in update_powers:
                functions[idx+power].append([idx])

        #convert to ANF specialized representation objects (easiest here):
        functions = [BooleanANF(ls) for ls in functions]

        # Mod by primitive polynomial P
        for idx in range(2*size-1, size-1, -1):
            for power in primitive_powers[:-1]:
                #XOR the ANFs
                functions[idx - size + power] ^= functions[idx]

        # convert to boolean functions and return the right subset
        self.fn_list = [fn.to_BooleanFunction() for fn in functions[:size]]

    @cached_property
    def minimal_polynomial(self) -> list[int]:
        """The minimal polynomial of the bit-0 output sequence.

        The state sequence of an MPR satisfies a linear recurrence over GF(2)
        whose characteristic polynomial divides the primitive polynomial P.
        For each individual bit, the observed sequence likewise satisfies a
        linear recurrence; its minimal polynomial is the lowest-degree such
        polynomial and equals the unique monic generator of the annihilator
        ideal of the sequence in GF(2)[x]. When the seed is in general
        position (as it is here, with seed value 2^n - 1) the minimal
        polynomial of bit 0 will typically equal P itself.

        Implementation: a length-2(n+1) prefix of the bit-0 sequence is
        observed and passed to Berlekamp-Massey, which is guaranteed to
        recover the minimal polynomial from any 2L consecutive symbols
        (where L is the linear complexity).

        :return: The coefficient list of the minimal polynomial, with index i
            holding the coefficient of x^i.
        :rtype: list[int]
        """
        #create a feedback register to determine minimal polynomial
        seq = []
        testReg = FeedbackRegister(2**self.size-1,self)
        for state in testReg.run((self.size+1)*2):
            seq.append(state[0])
        _, m = berlekamp_massey(seq)
        return [int(x) for x in m]