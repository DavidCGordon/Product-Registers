import numpy as np
import galois as gal

from itertools import product

# Rational Polynomial class for the entries of the matrix. Uses Galois GF(2) matrixrices.
class RationalPolynomial:
    @classmethod
    def unit(self):
        return RationalPolynomial([1],[1])
    @classmethod
    def zero(self):
        return RationalPolynomial([0],[1])
    
    @classmethod
    def from_int(self,value):
        return RationalPolynomial([value], [1])

    def __init__(self,n,d):
        if type(n) != gal.Poly:
            n = gal.Poly(n)
        if type(d) != gal.Poly:
            d = gal.Poly(d)

        self.n = n
        self.d = d

    def __add__(self,other):
        if type(other) != RationalPolynomial:
            raise ValueError(f"argument must be RationalPolynomial, not {type(other)}")
        out_n = self.d * other.n + self.n * other.d
        out_d = self.d * other.d
        return RationalPolynomial(out_n,out_d).simplify()

    def __mul__(self,other):
        if type(other) != RationalPolynomial:
            raise ValueError(f"argument must be RationalPolynomial, not {type(other)}")
        out_n = self.n * other.n
        out_d = self.d * other.d
        return RationalPolynomial(out_n,out_d).simplify()
    
    def __pow__(self,power):
        if type(power) != int:
            raise ValueError(f"power must be int, not {type(power)}")
        acc = RationalPolynomial.unit()
        for _ in range(power):
            acc *= self
        return acc

    def __truediv__(self,other):
        if type(other) != RationalPolynomial:
            raise ValueError(f"argument must be RationalPolynomial, not {type(other)}")
        out_n = self.n * other.d
        out_d = self.d * other.n
        return RationalPolynomial(out_n,out_d).simplify()

    def simplify(self):
        g = gal.gcd(self.n,self.d)
        return RationalPolynomial(self.n//g, self.d//g)
    
    def __str__(self):
        return "(" + str(self.n) + " / " + str(self.d) + ")"
    
    def __repr__(self):
        return str(self)

    def __eq__(self,other):
        other = self.from_int(other)
        return self.n == other.n and self.d == other.d


#technically would work for other fields?:

#creates an identity matrix given a field class:
def field_eye(field, size):
    return np.asarray(
        [field.unit() if i == j else field.zero() for i,j in product(range(size),repeat=2)],
        dtype = field
    ).reshape([size,size])

# Inverts a matrix, assuming it is square
def field_invert(field, matrix):
    #check to ensure matrix is square:
    if len(matrix.shape) != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"Matrix be square, not shape {matrix.shape}")
    size = matrix.shape[0]

    #append the identity, to be transformed into the inverse
    appended = field_eye(field, size)
    matrix = np.concatenate([matrix,appended], axis = 1)
    matrix.dtype = field

    #gaussian reduction:
    for pivot in range(size):
        # swap for a row with nonzero pivot
        offset = ((matrix[pivot:, pivot]) != 0).argmax()
        matrix[[pivot, pivot+offset]] = matrix[[pivot+offset, pivot]]

        # find the values we are zeroing out and create a matrix of changes
        pivot_value = matrix[pivot, pivot]
        row_multipliers = matrix[pivot+1:, pivot][:,np.newaxis] / pivot_value
        change = row_multipliers * matrix[pivot]

        # apply the changes to update matrix (+ and - are the same in GF(2))
        matrix[pivot+1:] += change

    # normalize rows by pivot
    for i in range(size):
        matrix[i] /= matrix[i][i]

    # backsubstitution/jordan reduction:
    for pivot in range(size-1,0,-1): 
        row_multipliers = matrix[:pivot, pivot][:,np.newaxis]
        change = row_multipliers * matrix[pivot]

        # apply the changes to update matrix (+ and - are the same in GF(2))
        matrix[:pivot] += change

    return (matrix[:, size:])