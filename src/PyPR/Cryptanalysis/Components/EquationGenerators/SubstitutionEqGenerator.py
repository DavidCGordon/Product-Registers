from typing import Any

from PyPR.BooleanLogic.FunctionInputs import VAR
from PyPR.BooleanLogic import BooleanFunction, BooleanANF

def SubstitutionEqGenerator(
    feedback_fn, 
    output_fn, 
    limit, 
):
    # Input handling:
    if type(output_fn) == list:
        return_list = True
        output_fn_list = output_fn
    else:
        return_list = False
        output_fn_list = [output_fn]

    bits = set.union(*(output_fn.idxs_used() for output_fn in output_fn_list))
    fns: list[Any] = [
        VAR(b) if b in bits else None 
        for b in range(feedback_fn.size)
    ]

    # Main Loop
    for t in range(limit+1):
        # yield equation (constant in ANF)
        if not return_list:
            yield (t, output_fn_list[0].compose(fns).translate_ANF(), 0)
        else:
            yield [
                (t, output_fn.compose(fns).translate_ANF(), 0)
                for output_fn in output_fn_list
            ]

        # don't do final update if not needed
        if t == (limit):
            break
        
        # update internal functions
        fns = [
            fns[b].compose(feedback_fn.fn_list).translate_ANF() 
            if b in bits else None
            for b in range(feedback_fn.size)
        ]
