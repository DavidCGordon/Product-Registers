from ProductRegisters.BooleanLogic import XOR, CONST
from ProductRegisters.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey, berlekamp_massey_iterator

import numpy as np
import galois as gal

from itertools import product
import re

# Rational Polynomial class for the entries of the matrix. Uses Galois GF(2) matrixrices.
class SequenceTransform:

    #several useful elements:
    @classmethod
    def one(self):
        return SequenceTransform([1],[1])
    @classmethod
    def zero(self):
        return SequenceTransform([0],[1])
    @classmethod
    def delay(self):
        return SequenceTransform([0,1],[1])
    
    # useful for broadcasting across arrays
    @classmethod
    def from_int(self,value):
        return SequenceTransform([value],[1])
    
    @classmethod
    def from_seq(self,seq):
        L,polynomial = berlekamp_massey(seq)
        arr = np.convolve(seq,polynomial)[:L+1] % 2
        return SequenceTransform(arr,polynomial)
    
    def __init__(self,n,d):
        if type(n) != gal.Poly:
            try:
                n = gal.Poly(n[::-1])
            except:
                print(n)
                raise ValueError(f"could not parse numerator input of type {type(n)} as a polynomial")
            
        if type(d) != gal.Poly:
            try:
                d = gal.Poly(d[::-1])
            except:
                print(d)
                raise ValueError(f"could not parse denominator input of type {type(d)} as a polynomial")

        self.n = n
        self.d = d


    def __add__(self,other):
        if type(other) != SequenceTransform:
            raise ValueError(f"argument must be SequenceTransform, not {type(other)}")
        out_n = self.d * other.n + self.n * other.d
        out_d = self.d * other.d
        return SequenceTransform(out_n,out_d).simplify()

    def __mul__(self,other):
        if type(other) != SequenceTransform:
            raise ValueError(f"argument must be SequenceTransform, not {type(other)}")
        out_n = self.n * other.n
        out_d = self.d * other.d
        return SequenceTransform(out_n,out_d).simplify()
    
    def __pow__(self,power):
        if type(power) != int:
            raise ValueError(f"power must be int, not {type(power)}")
        acc = SequenceTransform.one()
        for _ in range(power):
            acc *= self
        return acc

    def __truediv__(self,other):
        if type(other) != SequenceTransform:
            raise ValueError(f"argument must be SequenceTransform, not {type(other)}")
        out_n = self.n * other.d
        out_d = self.d * other.n
        return SequenceTransform(out_n,out_d).simplify()

    def simplify(self):
        g = gal.gcd(self.n,self.d)
        return SequenceTransform(self.n//g, self.d//g)
    
    # string formatting for z-transform is a bit of a pain :(
    def z_string(self):
        s = self.__str__().replace('D','z^(-1)')
        return re.sub(
            pattern = '\(-1\)\^(\d+)',
            repl = lambda x:  '(-' + x.group(1) + ')',
            string = s
        )
    
    def __str__(self):
        return (
            "(" + str(self.n).replace('x','D') + " / " + str(self.d).replace('x','D') + ")"
        )

    # so that it displays nicely in vectors
    def __repr__(self):
        return str(self)

    def __eq__(self,other):
        if type(other) != SequenceTransform:
            raise ValueError(f'Expected type SequenceTransform, not {type(other)}')
        return self.n == other.n and self.d == other.d
    
    def __copy__(self):
        return SequenceTransform(
            gal.Poly(self.n.coefficients()),
            gal.Poly(self.d.coefficients())
        )



# Methods for solving for the resolvent:

#creates an identity matrix given a field class:
def field_eye(field, size):
    entry_list = []
    for i,j in product(range(size),repeat=2):
        if i == j: entry_list.append(field.one())
        else: entry_list.append(field.zero())

    return np.asarray(entry_list,dtype=SequenceTransform).reshape([size,size])

# Gaussian Elimination matrix inversion:
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
        offset = ((matrix[pivot:, pivot]) != field.zero()).argmax()
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


# initialize useful variables:
def generate_resolvent_example(cmpr, use_z_convention = False):
    # Define variables
    D = SequenceTransform.delay()
    z = SequenceTransform.one() / SequenceTransform.delay()

    # collect information about CMPR:
    prev_state = cmpr._prev_state
    initial_state = cmpr._state
    REs = cmpr.fn.root_expressions()
    LC_bound = max([re.upper() for re in REs])

    # Initalize output dictionary 
    outputs = {
        'initial transforms': None,
        'sequence transforms': None,
        'chaining transforms': None,
        'scaled chaining': None,
        'scaled initial': None,
        'combined vector': None,
        'update matrices': cmpr.fn.update_matrices,
        'resolvent matrices': cmpr.fn.resolvent_matrices,
        'computed transforms':  np.array([None for i in range(cmpr.size)],dtype=SequenceTransform),
    }
    
    # fill in initial state transforms
    initial_vector = np.array(
        [SequenceTransform([cmpr[bit]],[1]) for bit in range(cmpr.size)],
        dtype=SequenceTransform
    )

    # iterate register to compute state/chaining sequences:
    register_values = [[0]*(2*LC_bound + 4) for bit in range(cmpr.size)]
    chaining_values = [[0]*(2*LC_bound + 4) for bit in range(cmpr.size)]
    chaining_fns = [XOR(CONST(0),*(bit_fn.args[1:])).compile() for bit_fn in cmpr.fn]

    cmpr.fn.compile()
    for t, state in enumerate(cmpr.run_compiled(2*LC_bound + 4)):
        for bit in range(cmpr.size):
            register_values[bit][t] = state[bit]
            chaining_values[bit][t] = chaining_fns[bit](state._state)

    register_vector = np.array([SequenceTransform.from_seq(seq) for seq in register_values])
    chaining_vector = np.array([SequenceTransform.from_seq(seq) for seq in chaining_values])

    outputs['sequence transforms'] = register_vector
    outputs['chaining transforms'] = chaining_vector
    outputs['initial transforms'] = initial_vector
    
    # compute everything using D-transform convention:
    # i.e. compute DC(D) + B[0]
    scaled_chaining = np.array([D]) * chaining_vector
    combined_vector = scaled_chaining + initial_vector
        
    outputs['scaled chaining'] = scaled_chaining
    outputs['combined vector'] = combined_vector

    # differences if using z-convention instead:
    if use_z_convention:
        # the initial, sequence, chaining, and computed transforms are all independent of notation:
        # using the 'z' convention used in the paper:
        
        # there is an additional delay = D = z^(-1) factor in the resolvents
        outputs['resolvent matrices'] = [np.array([D]) * matrix for matrix in cmpr.fn.resolvent_matrices]

        # and a 1/delay = 1/D = z factor in the combined vector.
        # when this distributed, it cancels the scaling on the chaining, and scales the initial instead
        del outputs['scaled chaining']
        outputs['scaled initial'] = np.array([z]) * initial_vector
        outputs['combined vector'] = np.array([z]) * combined_vector

        # in total, if resolvent is R and combined vector is v:
        # using this convention computes DR @ v/D instead of the equivalent R @ v
        # i.e. the effects of this modification simply cancel
    else:
        del outputs['scaled initial']

    # compute the desired product:
    for block_idx, block in enumerate(cmpr.fn.blocks):
        resolvent = outputs['resolvent matrices'][block_idx]
        input_vector = outputs['combined vector'][block]
        computed_vector = resolvent @ input_vector
        for bit in block:
            outputs['computed transforms'][bit] = computed_vector[bit-min(block)]

    # return register to its initial state:
    cmpr._prev_state = prev_state
    cmpr._state = initial_state
    return outputs

def pretty_print_example(output_dict, use_z_convention = False):
    # determine function used for printing:
    if use_z_convention: display = np.vectorize(lambda x: x.z_string())
    else: display = np.vectorize(lambda x: x.__str__())

    num_bits = len(output_dict['initial transforms'])
    num_blocks = len(output_dict['update matrices'])

    for output_name, value in output_dict.items():
        print(f"{output_name}:")
        if output_name == 'update matrices':
            for block in range(num_blocks):
                print(f"Component {block}:\n{(value[block])}")
        elif output_name == 'resolvent matrices':
            for block in range(num_blocks):
                print(f"Component {block}:\n{display(value[block])}")
        else:
            for bit in range(num_bits):
                print(f"Bit {bit}: {display(value[bit])}")
        print('\n')