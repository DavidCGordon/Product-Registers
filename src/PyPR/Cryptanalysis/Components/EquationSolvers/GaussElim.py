import numpy as np
import numba

from PyPR.Cryptanalysis.Components.EquationStores.EqStore import EqStore
from PyPR.Cryptanalysis.Components.EquationStores.LUEqStore import LUEqStore
from PyPR.Cryptanalysis.Components.EquationStores.LUDynamicEqStore import LUDynamicEqStore
from PyPR.Cryptanalysis.Components.EquationStores.DynamicEqStore import DynamicEqStore

u8 = numba.types.uint8
@numba.njit(numba.types.Tuple((u8[:,:],u8[:]))(u8[:,:]))
def reduce_matrix(matrix):
    # returns the rref and the free vars
    rows, cols = matrix.shape
    p_row, p_col = 0, 0
    free_vars = np.zeros((cols,),dtype='uint8')
    
    while p_row < rows and p_col < cols:
        # Find the pivot element/swap rows
        for i in range(p_row + 1, rows):
            if matrix[i,p_col] > matrix[p_row,p_col]:
                matrix[np.array([p_row, i])] = matrix[np.array([i, p_row])]
                break

        # Identify free vars
        if matrix[p_row,p_col] == 0:
            free_vars[p_col] = 1
            p_col += 1
            continue
        
        # No normalization needed due to 0/1 only

        # Eliminate other rows
        for i in range(rows):
            if i != p_row and matrix[i,p_col]:
                matrix[i] ^= matrix[p_row]

        p_row += 1
        p_col += 1

    # make sure the last columns are counted as free:
    for i in range(p_col,cols):
        free_vars[i] = 1
 
    return matrix[:p_row], free_vars



def solve(equation_store, constants = None):
    if type(equation_store) in [EqStore, DynamicEqStore]:
        matrix = equation_store.eqs
    else:
        raise ValueError
    
    rref, free_vars = reduce_matrix(matrix)
