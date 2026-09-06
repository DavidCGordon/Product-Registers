from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey

import galois as gl
import numpy as np
import re

# Rational Polynomial class for the entries of the matrix. Uses Galois GF(2) matrices.
class BooleanGF:

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
        if type(other) != BooleanGF:
            raise ValueError(f'Expected type BooleanGF, not {type(other)}')
        return self.num == other.num and self.den == other.den
    
    def __copy__(self):
        return BooleanGF(
            gl.Poly(self.num.coefficients()),
            gl.Poly(self.den.coefficients())
        )
    
    @classmethod
    def from_seq(cls,seq):
        L,polynomial = berlekamp_massey(seq)
        arr = np.convolve(seq,polynomial)[:L+1] % 2
        return BooleanGF(arr,polynomial)
    
    def to_seq(self):
        pass