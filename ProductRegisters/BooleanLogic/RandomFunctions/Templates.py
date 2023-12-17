from ProductRegisters.BooleanLogic import *

import random

def uniform(n):
    return [1/n for i in range(n)]

def geometric(r,n,c = None):
    if c: return [(c*(r**i)) for i in range(n)]  
    else: return [(1-r) * (r**i) / (1-r**(n)) for i in range(n)]

def constant(c,n):
    return [c for i in range(n)]


def old_ANF_chaining_template(max_vars = 4, max_terms = 4, individual_terms = True, added_constants = False):
    def add_constants(cmpr,fn, allowed_blocks):
        if random.random() > .5:
            fn.arg_limit += 1
            fn.add_arguments(CONST(1))
        return (True,fn)

    def default_modifier(cmpr,fn,allowed_blocks):
        return (True, fn)

    def template(cmpr):
        kwargs = {}
        kwargs["max_depth"] = 2
        kwargs["current_depth"] = 0 

        n = cmpr.num_components
        kwargs["block_probabilities"] = [len(block)/cmpr.size for block in cmpr.blocks]
        kwargs["input_densities"] = constant(1,n)
        kwargs["input_minimums"] = constant(1,n)

        if individual_terms:
            depth_0 = {XOR(arg_limit = i): 1/(max_terms) for i in range(1,max_terms+1)}
            depth_1 = {AND(arg_limit = i): 1/(max_vars) for i in range(1,max_vars+1)}
        else:
            depth_0 = {XOR(arg_limit = i): 1/(max_terms-1) for i in range(2,max_terms+1)}
            depth_1 = {AND(arg_limit = i): 1/(max_vars-1) for i in range(2,max_vars+1)}
            
        distributions = [[depth_0,depth_1]]
        kwargs["component_distributions"] = distributions
        kwargs["depth_mode"] = "discrete"

        if added_constants:
            kwargs["modifier"] = add_constants
        else:
            kwargs["modifier"] = default_modifier

        return kwargs
    return template



def and_or_mixture_template(max_vars = 4, max_terms = 4):
    def template(cmpr):
        kwargs = {}
        kwargs["max_depth"] = 3
        kwargs["current_depth"] = 0 

        n = cmpr.num_components
        square_sum = sum(len(block) for block in cmpr.blocks)
        kwargs["block_probabilities"] = [len(block)/square_sum for block in cmpr.blocks]
        kwargs["input_densities"] = constant(.5,n)
        kwargs["input_minimums"] = constant(1,n)

        depth_0 = {XOR(arg_limit = i): 1/(max_terms) for i in range(1,max_terms+1)}
        depth_1 = {
            **{AND(arg_limit = i): .2/(max_vars) for i in range(1,max_vars+1)},
            **{OR(arg_limit = i): .7/(max_vars) for i in range(1,max_vars+1)},
            **{VAR(None): .05 ,CONST(1): .05}
        }
        depth_2 = {
            **{AND(arg_limit = i): .7/(max_vars) for i in range(1,max_vars+1)},
            **{OR(arg_limit = i): .2/(max_vars) for i in range(1,max_vars+1)},
            **{VAR(None): .05 ,CONST(1): .05}
        }
            
        distributions = [[depth_0,depth_1,depth_2]]
        kwargs["component_distributions"] = distributions
        kwargs["depth_mode"] = "discrete"

        return kwargs
    return template