from ProductRegisters.BooleanLogic import XOR,AND,OR,NOR,NAND,NOT,CONST,VAR
from numpy.random import choice
from random import sample, random
from math import ceil


def leaves(component, depth_inc = 1):
    if hasattr(component, "args"):
        for arg in component.args:
            for leaf, depth_inc in leaves(arg, depth_inc = depth_inc + 1):
                yield leaf, depth_inc 

        if len(component.args) < component.arg_limit:
            yield component, depth_inc


#default components selected for ANF style gates:
default_components = [ AND(arg_limit = 4), AND(arg_limit = 3), AND(arg_limit = 2),
                       XOR(arg_limit = 4), XOR(arg_limit = 3), XOR(arg_limit = 2),
                       CONST(1), VAR(None),
]


def random_function(
        allowed_blocks,
        current_depth = None,
        max_depth = None,
        block_probabilities = None,
        component_distributions = None,
        depth_mode = None
    ):

        # handle the bottom layer case
        if current_depth >= max_depth: 
            return VAR(0)

        # or choose a component
        else:
            # get components and probabilities from the appropriate dict
            distribution_idx = min(len(component_distributions)-1, current_depth)
            distribution_list = list(component_distributions[distribution_idx].items())
            allowed_components = [x[0] for x in distribution_list]
            probabilities = [x[1] for x in distribution_list]

            # create the component
            component_choice = choice(
                a = allowed_components,
                p = probabilities
            ).__copy__()
        

        # Handle VAR and CONST components
        if type(component_choice) == CONST:
            return component_choice
        elif type(component_choice) == VAR:
            # return VAR(choice(allowed_bits))
            return component_choice


        # fill out any dangling leaves with recursive calls
        for leaf, depth_inc in leaves(component_choice):
            
            num_children = leaf.arg_limit - len(leaf.args)
            
            # update depth for recursive calls.
            if depth_mode == "discrete":
                new_depth = current_depth + 1
            elif depth_mode == "granular":
                new_depth = current_depth + depth_inc
            else:
                raise ValueError('depth_mode must be either "discrete" or "granular".')

            
            #information used to ensure children are selected w/o repetition
            children = []
            vars_used = set()
            blocks_used = set()

            # recursive calls to generate the additional children for this leaf
            for _ in range(num_children):
                new_child = random_function(
                    allowed_blocks = allowed_blocks,
                    current_depth = new_depth,
                    max_depth = max_depth,
                    block_probabilities = block_probabilities,
                    component_distributions = component_distributions,
                    depth_mode = depth_mode
                )

                if type(new_child) == VAR:
                    # establish which blocks/probabilites are allowed:
                    blocks_available = []
                    probs_available = []

                    for block_idx in range(len(allowed_blocks)):
                        if not (block_idx in blocks_used):
                            blocks_available.append(block_idx)
                            probs_available.append(block_probabilities[block_idx])
                    
                    probs_available = [x / sum(probs_available) for x in probs_available]

                    # if there are is no space_left (in theory this should never happen):
                    if not blocks_available: 
                        #children.append(VAR(choice(list(set().union(*allowed_blocks)))))
                        leaf.arg_limit -= 1
                        continue 
                    else:
                        block_choice = int(choice(blocks_available, p = probs_available))

                    # choose bit
                    bits_available = list(set(allowed_blocks[block_choice]) - vars_used)
                    bit_choice = int(choice(bits_available))

                    # update used variables/blocks sets
                    vars_used.add(bit_choice)
                    if all(x in vars_used for x in allowed_blocks[block_choice]):
                        blocks_used.add(block_choice)
                    new_child = VAR(bit_choice)
                else:
                    # TODO: eliminate duplicate branches
                    # might require an efficient equals function, which may be computationally difficult.
                    pass

                # append new child:
                children.append(new_child)
            leaf.add_arguments(*children)
        return component_choice

