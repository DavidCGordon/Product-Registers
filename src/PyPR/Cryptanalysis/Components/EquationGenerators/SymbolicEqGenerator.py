from PyPR.BooleanLogic import BooleanANF

def SymbolicEqGenerator(
    feedback_fn, 
    output_fn, 
    limit, 
    initialization=None
):
    # Input handling:
    if type(output_fn) == list:
        return_list = True
        output_fn_list = output_fn
    else:
        return_list = False
        output_fn_list = [output_fn]

    if initialization == None:
        fns = [BooleanANF([[b]]) for b in range(feedback_fn.size)]
    else:
        fns = [BooleanANF.from_BooleanFunction(f) for f in initialization]

    # set up the functions to evaluate
    eval_list = [
        f.anf_optimize().remap_constants([
            (0, BooleanANF()),
            (1, BooleanANF([True]))
        ]) for f in feedback_fn.fn_list
    ]

    output_fn_list = [
        f.anf_optimize().remap_constants([
            (0, BooleanANF()),
            (1, BooleanANF([True]))
        ]) for f in output_fn_list
    ]

    # main loop
    for t in range(limit+1):
        # yield equation (constant in ANF
        if not return_list:
            yield (t, output_fn_list[0].eval_ANF(fns).to_BooleanFunction(), 0)
        else:
            yield [
                (t, output_fn.eval_ANF(fns).to_BooleanFunction(), 0)
                for output_fn in output_fn_list
            ]

        # don't do final update if not needed
        if t == (limit):
            break
        
        # update internal functions
        fns = [eval_list[b].eval_ANF(fns) for b in range(feedback_fn.size)]
