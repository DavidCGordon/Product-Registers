import numpy as np
import numba

from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore

u8 = numba.types.uint8
@numba.njit(u8[:](u8[:,:],u8[:,:],u8[:]))
def lu_solve(L,U,z):

    # backsolve L: z  =>  L^{-1}z
    for i in range(len(z)-1):
        for j in range(i+1,len(z)):
            z[j] ^= L[j,i] * z[i]

    # backsolve U: L^{-1}z => U^{-1}L^{-1}z
    # equivalently: (LU)^{-1}(z + c)
    for i in range(len(z)-1,0,-1):
        for j in range(i):
            z[j] ^= U[j,i] * z[i]

    return z

# If eq store is inconsistent, you might have additional constants:
def solve(equation_store, additional_constants = None):
    """
    Solve with specific constants.

    things are zeroed by default,

    :param equation_store: _description_
    :type equation_store: _type_
    :param additional_constants: _description_, defaults to None
    :type additional_constants: _type_, optional
    :raises ValueError: _description_
    :raises ValueError: _description_
    :raises ValueError: _description_
    :return: _description_
    :rtype: _type_
    """
    if type(equation_store) != LUEqStore:
        raise ValueError("LU Solver only takes LU equation stores")
    
    if (
        additional_constants is not None and 
        len(additional_constants) != (equation_store.num_vars)
    ):
        raise ValueError("Bad Length")
    

    if equation_store.consistent:
        if additional_constants is None:
            constant_vector = np.zeros([equation_store.num_vars], dtype=np.uint8)
        else:
            constant_vector = additional_constants.copy()
        constant_vector[equation_store.comb_to_idx[tuple()]] = 1

    elif not equation_store.consistent:
        if additional_constants is None:
            constant_vector = np.zeros([equation_store.num_vars], dtype=np.uint8)
        else:
            constant_vector = additional_constants.copy()

    n = equation_store.num_vars
    sol = lu_solve(
        equation_store.lower_matrix[:n,:n],
        equation_store.upper_matrix[:n,:n],
        constant_vector,
    )

    return sol

def solve_with_guesses(equation_store):
    pass # TODO: have a version with the guessing / pruning logic here :) 
         # generator!!
         # do