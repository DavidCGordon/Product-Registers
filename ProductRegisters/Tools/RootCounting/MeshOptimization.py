# Linear Complexity and Monomial estimation
from ProductRegisters.Tools.RootCounting.MonomialProfile import TermSet,MonomialProfile
from ProductRegisters.Tools.RootCounting.JordanSet import JordanSet
from ProductRegisters.Tools.RootCounting.RootExpression import RootExpression

# Other libs
import numpy as np
from numba import njit
from numba.types import int32 as i32

# use the iterator and convert output to a RE
def compute_single_mesh(sizes,degrees,locked_list=None):
    deg_arr = np.array(degrees, dtype='int32')
    size_arr = np.array(sizes, dtype='int32')
    if not locked_list:
        locked_arr = np.ones_like(size_arr)
    else:
        locked_arr = np.array(locked_list, dtype = 'int32')

    assert (len(sizes) == len(locked_list))

    root_table = {}
    for output in _re_mesh_iterator(size_arr,deg_arr,locked_arr):
        # convert tuple to jset
        field = []
        counts = []
        for size, count in zip(size_arr,output):
            if count != 0:
                field.append(size)
                counts.append(count)
        field_tuple = tuple(sorted(field))
        js = JordanSet({k:min(k,v) for k,v, in zip(field,counts)}, set([1]))

        # add to the table
        if field_tuple in root_table:
            root_table[field_tuple].add(js)
        else: 
            root_table[field_tuple] = set([js])

    return RootExpression(root_table)



# return arrays representing all the JordanSets
@njit(i32[:](i32[:],i32[:],i32[:]))
def _re_mesh_iterator(sizes, degrees, locked_list):
    sizes *= locked_list

    # create the values matrix:
    values = degrees.copy()
    for i in range(1,len(values)):
        values[i] = degrees[i] * values[i-1]

    arr = np.zeros((len(degrees),len(degrees)),dtype='int32')
    arr[0][len(degrees)-1] = 1

    depth = 0

    yield arr[depth]
    while depth >= 0:
        block = (len(degrees)-1)-depth

        # return condition
        max_depth_reached = (depth == len(degrees)-1)
        no_more_trades = (arr[depth][block] == 0)
        #end_saturated = np.all((arr[depth] > sizes)[:block])
        if max_depth_reached or no_more_trades: # or end_saturated:
            depth -= 1
            continue
        
        # update using the given trades
        num_traded = max(1, arr[depth][block] - sizes[block])
        arr[depth][block] -= num_traded
        arr[depth][block-1] += num_traded * degrees[block]


        # yield before moving down tree (pre-order)
        saturated_mask = (arr[depth] >= sizes)
        partially_filled_mask = ((0 < arr[depth]) & (arr[depth] < sizes))

        # ensure uniqueness in the assignments to the saturated field by only
        # yielding assignments where all extra is in the leftmost saturated field
        leftmost_assignment = True
        first_idx_already_found = False
        for i in range(len(degrees)):
            if saturated_mask[i]:
                if first_idx_already_found:
                    if (arr[depth][i] > sizes[i]):
                        leftmost_assignment = False
                        break
                else:
                    first_idx_already_found = True
        
        # check additionally, that any partially filled field can't be incremented further. 
        # if it can, then that assignment is strictly better, so don't yield this one. 
        fillable_space = False
        if leftmost_assignment and np.any(partially_filled_mask):
            overshoot = np.maximum(0,arr[depth] - sizes)
            extra_potential = np.sum(overshoot * values) / values
            potential_increment = extra_potential * partially_filled_mask.astype('int32')
            fillable_space = np.any(potential_increment >= 1)
        

        if leftmost_assignment and not fillable_space:
            yield arr[depth]


        # recurse down tree
        arr[depth + 1] = arr[depth]
        depth += 1

    raise StopIteration()