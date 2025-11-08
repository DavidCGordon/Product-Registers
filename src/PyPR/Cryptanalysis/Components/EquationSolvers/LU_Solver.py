import numpy as np
import numba

from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore
from PyPR.Cryptanalysis.Components.EquationStores.LUDynamicEqStore import LUDynamicEqStore
from PyPR.Cryptanalysis.Components.EquationStores.DynamicEqStore import DynamicEqStore

u8 = numba.types.uint8
@numba.njit(u8[:](u8[:,:],u8[:,:],u8[:],u8[:]))
def lu_solve(L,U,c,z):
    z = z.copy()

    # backsolve L: z  =>  L^{-1}z
    for i in range(len(z)-1):
        for j in range(i+1,len(z)):
            z[j] ^= L[j,i] * z[i]

    # add in constants: L^{-1}z  =>  L^{-1}(z + c)
    # (constants are stored as L^{-1}c)
    for i in range(len(z)):
        z[i] ^= c[i]

    # backsolve U: L^{-1}(z + c) => U^{-1}L^{-1}(z + c)
    # equivalently: (LU)^{-1}(z + c)
    for i in range(len(z)-1,0,-1):
        for j in range(i):
            z[j] ^= U[j,i] * z[i]
            
    return z

# Incomplete solver -> additional constants
def solve(equation_store, additional_constants = None):
    if additional_constants is None:
        additional_constants = np.zeros_like(equation_store.constants)
    if type(equation_store) in [LUEqStore, LUDynamicEqStore]:
        return lu_solve(
            equation_store.lower_matrix,
            equation_store.upper_matrix,
            equation_store.constants,
            additional_constants
        )
    else:
        raise ValueError("LU Solver only takes LU equation stores")
    