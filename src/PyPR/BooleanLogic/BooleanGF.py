import re
from typing import Any, Self

import galois as gl
import numpy as np

from PyPR.JSON_Serialization import Serializable

from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey


# Rational Polynomial class for the entries of the matrix. Uses Galois GF(2) matrices.
class BooleanGF(Serializable):

    #several useful elements:
    @classmethod
    def one(cls):
        return BooleanGF([1],[1])
    @classmethod
    def zero(cls):
        return BooleanGF([0],[1])
    @classmethod
    def delay(cls):
        return BooleanGF([0,1],[1])

    # useful for broadcasting across arrays
    @classmethod
    def from_int(cls,value):
        return BooleanGF([value],[1])

    def __init__(self,numerator,denominator):
        # Sanitize Numerator:
        if type(numerator) != gl.Poly:
            try:
                numerator = gl.Poly(numerator[::-1])
            except:
                raise ValueError(f"could not parse numerator input of type {type(numerator)} as a polynomial")

        # Sanitize Denomintaor:
        if type(denominator) != gl.Poly:
            try:
                denominator = gl.Poly(denominator[::-1])
            except:
                print(denominator)
                raise ValueError(f"could not parse denominator input of type {type(denominator)} as a polynomial")

        self.num = numerator
        self.den = denominator

    def simplify(self):
        g = gl.gcd(self.num,self.den)
        return BooleanGF(self.num//g, self.den//g)

    def __add__(self,other):
        if type(other) != BooleanGF:
            raise ValueError(f"argument must be BooleanGF, not {type(other)}")
        out_n = self.den * other.num + self.num * other.den
        out_d = self.den * other.den
        return BooleanGF(out_n,out_d).simplify()

    def __mul__(self,other):
        if type(other) != BooleanGF:
            raise ValueError(f"argument must be BooleanGF, not {type(other)}")
        out_n = self.num * other.num
        out_d = self.den * other.den
        return BooleanGF(out_n,out_d).simplify()

    def __pow__(self,power):
        if type(power) != int:
            raise ValueError(f"power must be int, not {type(power)}")
        acc = BooleanGF.one()
        for _ in range(power):
            acc *= self
        return acc

    def __truediv__(self,other):
        if type(other) != BooleanGF:
            raise ValueError(f"argument must be BooleanGF, not {type(other)}")
        out_n = self.num * other.den
        out_d = self.den * other.num
        return BooleanGF(out_n,out_d).simplify()

    # string formatting for z-transform is a bit of a pain :(
    def _z__str__(self):
        s = self.__str__().replace('D','z^(-1)')
        return re.sub(
            pattern = "\\(-1\\)\\^(\\d+)",
            repl = lambda x:  '(-' + x.group(1) + ')',
            string = s
        )

    def __str__(self):
        return (
            "(" + str(self.num).replace('x','D') + " / " + str(self.den).replace('x','D') + ")"
        )

    # so that it displays nicely in vectors
    def __repr__(self):
        return str(self)

    def __eq__(self,other):
        # NotImplemented rather than an error for a foreign operand, as the data
        # model requires: Python compares dict keys and set members of any type
        # with ==, and serialization keys its id map by object
        if not isinstance(other, BooleanGF):
            return NotImplemented
        return self.num == other.num and self.den == other.den

    def __hash__(self):
        # equal exactly when both coefficient lists are, so hash those
        return hash((
            tuple(self.num.coefficients().tolist()),
            tuple(self.den.coefficients().tolist()),
        ))

    def __copy__(self):
        return BooleanGF(
            gl.Poly(self.num.coefficients()),
            gl.Poly(self.den.coefficients())
        )

    # Serialization (the id hook is the Serializable default: a BooleanGF refers
    # to no other serializable object, and it is never modified in place)
    def _generate_JSON_entry(self,
        ids: dict[Any, int]
    ) -> dict[str, Any]:
        """Write the numerator and denominator as coefficient lists.

        The lists use the constructor's own order -- index `i` holds the
        coefficient of `D^i` -- so parsing is a constructor call. The field is
        not stored: it must be GF(2), which is what the constructor builds.

        :param ids: The map from objects to ids; unused.
        :type ids: dict[Any, int]
        :raises ValueError: If either polynomial is over a field other than GF(2).
        :return: The fraction's data.
        :rtype: dict[str, Any]
        """
        for poly in (self.num, self.den):
            if poly.field.order != 2:
                raise ValueError(
                    f"only BooleanGF over GF(2) can be written to JSON, not over {poly.field.name}"
                )
        return {
            "numerator": self.num.coefficients()[::-1].tolist(),
            "denominator": self.den.coefficients()[::-1].tolist(),
        }

    @classmethod
    def _parse_JSON_entry(cls,
        object_data: dict[str, Any],
        parsed_objects: list[Any]
    ) -> Self:
        """Rebuild a BooleanGF from the coefficient lists `_generate_JSON_entry` wrote.

        :param object_data: The data written for this fraction.
        :type object_data: dict[str, Any]
        :param parsed_objects: The objects rebuilt so far; unused.
        :type parsed_objects: list[Any]
        :return: The rebuilt fraction.
        :rtype: Self
        """
        return cls(object_data["numerator"], object_data["denominator"])

    @classmethod
    def from_seq(cls,seq):
        L,polynomial = berlekamp_massey(seq)
        arr = np.convolve(seq,polynomial)[:L+1] % 2
        return BooleanGF(arr,polynomial)

    def to_seq(self):
        pass
