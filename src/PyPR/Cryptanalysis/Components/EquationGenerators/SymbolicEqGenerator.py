from typing import Iterator
from PyPR.BooleanLogic import BooleanANF, BooleanFunction
from PyPR.FeedbackFunctions import FeedbackFunction

def SymbolicEqGenerator(
    feedback_fn: FeedbackFunction, 
    output_fn: BooleanFunction | list[BooleanFunction], 
    limit, 
    initialization=None
) -> Iterator[
    BooleanFunction |
    list[BooleanFunction]
]:
    """Generates the Equations for a given feedback function and output function. Generally this is
    slightly slower than the symbolic equation generator, and dramatically slower than the cube
    equation generator.

    Like all equation generators, the return is an iterator which iterates through tuples of the
    form `(time, equation, extra_constant)`. Time is an integer which denotes which clock cycle the
    equation represents. The equation can be different types depending on the equation generator,
    but for the SymbolicEqGenerator, equations are generated as BooleanFunctions. The extra constant 
    is always 0, as the constant is included in the BooleanFunction as a `CONST` node. 
    
    This function also allows you to pass in a list of output functions. If you do this,
    the, output which is yielded on each cycle will be a list of tuples, with each corresponding
    to one of the output functions passed in. These tuples each have the same (time, equation, 
    extra_constant) format as if only one output function is passed.\n

    ## How it works:
    On the `n`th iteration, the stored equations represent the equations at time `t=n`, in terms of
    the variables at time `t=0`. We can then simply compose the feedback_function equations with
    the current equations to get the equations for time `t=n+1`. `eval_ANF` on the composed functions
    simply evaluates the ANF operations symbolically.
    
    ## Comparison to Substitution Equation Generator:
    In implementation, these two generators are almost identical. The main line which is different,
    which replaces the functions, simply composes in a different direction:
     - **Symbolic:** `current_equations = feedback_functions.compose(current_equations)`
     - **Substitution:** `current_equations = current_equations.compose(feedback_functions)`

    However, there's actually a surprising number of consequences to this:
    - The symbolic equation generator allows you to give initial symbolic expressions (useful for 
        describing on of the starting state variables in terms of initialization, or replacing known
        variables with constants), and works correctly when manipulating symbols. Because the
        substition generator assumes the `t=0` and `t=1` equations in order to work, it doesnt have 
        this ability. This can drastically simplify the equations (compared to generating them in full 
        complexity and then simplifying) and if you need/can take advantage of this functionality, 
        you should use the symbolic generator.
     - In the substitution generator, we only need the `t=n` equation for a bit to compute the `t=n+1`
        equation for the same bit. In contrast, the symbolic equation generator needs expressions for
        every bit. This means that if only a small portion of the bits are needed to compute the output
        function (which is often the case), the substituion generator can do less work by only computing
        equations for those bits, making it generally faster
     - Lastly, it changes the size distributions of the ANF multiplications involved,
        which can affect performance. This is fairly subtle and complicated, so it's discussed 
        in more depth in the performance section

    ## Performance:
    The main speed gains come from a couple of things:
    1. There is huge amount of reuse, both in reusing evalutions of the states in different cubes, and
        reusing older coefficients to cover large amounts of the sum. 
    2. The summation indices are the same everytime; after precomputation, this looks like a very tight loop
        of summing values in an array. The storage of everything in dense arrays also helps with cache locality
    3. because the algorithm is so computer-friendly (especially relative to handling symbolics), we can
        numba compile it to speed it up even further. 
    altogether this is probably 2-3 orders of magnitude faster than the other methods

    :param feedback_fn: The function used to update the state on each clock cycle
    :type feedback_fn: FeedbackFunction
    :param output_fn: The function used to derive output from the state
    :type output_fn: BooleanFunction
    :param limit: The number of cycles to generate output for
    :type limit: int
    :param verbose: whether or not to print output, defaults to False
    :type verbose: bool, optional
    :param _print_depth: The indentation level to print at (advise not to touch this, it's 
        mostly internal to make printing prettier, and doesnt change much), defaults to 0
    :type _print_depth: int, optional
    :return: An Iterator which yields (time,equation,extra_const) tuples. If a list of output 
        functions are passed as input, the iterator will yield a list of such Tuples on each iteration.
        Otherwise, only one tuple will be yielded each iteration. 
    :rtype: Iterator[
        tuple[int, BooleanFunction, int] | 
        list[tuple[int, BooleanFunction, int]] 
    ]
    """
    # set flags to match outputs shape to input shape:
    if type(output_fn) == list:
        return_list = True
        output_fn_list = output_fn
    elif isinstance(output_fn,BooleanFunction):
        return_list = False
        output_fn_list = [output_fn]
    else:
        raise TypeError(
            f"output_fn must be a BooleanFunction or list of Boolean functions. " +
            f"Got {type(output_fn)} instead."
        )
    
    if initialization == None:
        fns = [BooleanANF([[b]]) for b in range(feedback_fn.size)]
    else:
        fns = [BooleanANF.from_BooleanFunction(f) for f in initialization]

    # set up the functions to evaluate
    eval_list = [
        #f.anf_optimize().remap_constants([
        f.remap_constants([
            (0, BooleanANF()),
            (1, BooleanANF([True]))
        ]) for f in feedback_fn.fn_list
    ]

    output_fn_list = [
        f.remap_constants([
            (0, BooleanANF()),
            (1, BooleanANF([True]))
        ]) for f in output_fn_list
    ]

    # main loop
    for t in range(limit+1):
        # yield equation (constant in ANF
        if not return_list:
            yield output_fn_list[0].eval_ANF(fns).to_BooleanFunction()
        else:
            yield [output_fn.eval_ANF(fns).to_BooleanFunction() for output_fn in output_fn_list]

        # don't do final update if not needed
        if t == (limit):
            break
        
        # update internal functions
        fns = [eval_list[b].eval_ANF(fns) for b in range(feedback_fn.size)]
