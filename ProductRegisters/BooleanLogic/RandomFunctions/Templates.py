from ProductRegisters.BooleanLogic import XOR,AND,CONST

import random

def uniform(n):
    return [1/n for i in range(n)]

def geometric(r,n,c = None):
    if c: return [(c*(r**i)) for i in range(n)]  
    else: return [(1-r) * (r**i) / (1-r**(n)) for i in range(n)]

def constant(c,n):
    return [c for i in range(n)]

# true values/individual terms
def old_ANF_chaining_template(max_terms = 4, max_ands = 5):
    def old_ANF_modifier(fn, allowed_blocks):
        if random.random() > .5:
            fn.arg_limit += 1
            fn.add_arguments(CONST(1))
        return fn

    def template(cmpr):
        kwargs = {}
        kwargs["max_depth"] = 2
        kwargs["current_depth"] = 0 

        n = cmpr.num_components
        kwargs["block_probabilities"] = [len(block)/cmpr.size for block in cmpr.blocks]
        kwargs["input_densities"] = constant(1,n)
        kwargs["input_minimums"] = constant(1,n)

        depth_0 = {XOR(arg_limit = i): 1/(max_ands) for i in range(1,max_ands+1)}
        depth_1 = {AND(arg_limit = i): 1/(max_terms) for i in range(1,max_terms+1)}
        
        distributions = [[depth_0,depth_1]]
        kwargs["component_distributions"] = distributions
        kwargs["depth_mode"] = "discrete"

        kwargs["modifier"] = old_ANF_modifier

        return kwargs
    return template

# No individual terms
def old_ANF_chaining_template2(max_terms = 4, max_ands = 5):
    def old_ANF_modifier(fn, allowed_blocks):
        if random.random() > .5:
            fn.arg_limit += 1
            fn.add_arguments(CONST(1))
        return fn

    def template(cmpr):
        kwargs = {}
        kwargs["max_depth"] = 2
        kwargs["current_depth"] = 0 

        n = cmpr.num_components
        kwargs["block_probabilities"] = [len(block)/cmpr.size for block in cmpr.blocks]
        kwargs["input_densities"] = constant(1,n)
        kwargs["input_minimums"] = constant(1,n)

        depth_0 = {XOR(arg_limit = i): 1/(max_terms-1) for i in range(2,max_terms+1)}
        depth_1 = {AND(arg_limit = i): 1/(max_terms-1) for i in range(2,max_terms+1)}
        
        distributions = [[depth_0,depth_1]]
        kwargs["component_distributions"] = distributions
        kwargs["depth_mode"] = "discrete"

        kwargs["modifier"] = old_ANF_modifier

        return kwargs
    return template

# no True values
def old_ANF_chaining_template3(max_terms = 4, max_ands = 5):
    def template(cmpr):
        kwargs = {}
        kwargs["max_depth"] = 2
        kwargs["current_depth"] = 0 

        n = cmpr.num_components
        kwargs["block_probabilities"] = [len(block)/cmpr.size for block in cmpr.blocks]
        kwargs["input_densities"] = constant(1,n)
        kwargs["input_minimums"] = constant(1,n)

        depth_0 = {XOR(arg_limit = i): 1/(max_terms) for i in range(1,max_terms+1)}
        depth_1 = {AND(arg_limit = i): 1/(max_terms) for i in range(1,max_terms+1)}
        
        distributions = [[depth_0,depth_1]]
        kwargs["component_distributions"] = distributions
        kwargs["depth_mode"] = "discrete"
        return kwargs
    return template

# No individual terms, no True values
def old_ANF_chaining_template4(max_terms = 4, max_ands = 5):
    def template(cmpr):
        kwargs = {}
        kwargs["max_depth"] = 2
        kwargs["current_depth"] = 0 

        n = cmpr.num_components
        kwargs["block_probabilities"] = [len(block)/cmpr.size for block in cmpr.blocks]
        kwargs["input_densities"] = constant(1,n)
        kwargs["input_minimums"] = constant(1,n)

        depth_0 = {XOR(arg_limit = i): 1/(max_terms-1) for i in range(2,max_terms+1)}
        depth_1 = {AND(arg_limit = i): 1/(max_terms-1) for i in range(2,max_terms+1)}
        
        distributions = [[depth_0,depth_1]]
        kwargs["component_distributions"] = distributions
        kwargs["depth_mode"] = "discrete"
        return kwargs
    return template


def better_ANF_chaining_template(max_terms = 4, max_ands = 4):
    def template(cmpr):
        kwargs = {}
        kwargs["max_depth"] = 2
        kwargs["current_depth"] = 0 

        n = cmpr.num_components
        kwargs["block_probabilities"] = geometric(1.3, n)
        kwargs["input_densities"] = geometric(.9, n, c = 1)
        kwargs["input_minimums"] = constant(1, n)

        depth_0 = {XOR(arg_limit = i): 1/(max_terms-1) for i in range(2,max_terms+1)}
        depth_1 = {AND(arg_limit = i): 1/(max_terms-1) for i in range(2,max_terms+1)}
        distributions = [[depth_0,depth_1]]
        kwargs["component_distributions"] = distributions
        kwargs["depth_mode"] = "discrete"

        return kwargs
    return template