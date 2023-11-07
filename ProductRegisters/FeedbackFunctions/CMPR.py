# Simulation Boilerplate
from ProductRegisters.BooleanLogic import ANF_spec_repr, XOR, AND, CONST, VAR
from ProductRegisters.FeedbackFunctions import FeedbackFunction
from ProductRegisters.FeedbackFunctions import MPR

# Chaining Generation/Handling:
from ProductRegisters.BooleanLogic.RandomFunctions import random_function

# Linear Complexity estimation
from ProductRegisters.Tools.RootCounting.JordanSet import JordanSet
from ProductRegisters.Tools.RootCounting.RootExpression import RootExpression

# Other analysis
import ProductRegisters.Tools.ResolventSolving as ResolventSolving
from ProductRegisters.Tools.MersenneTools import expected_period, expected_period_ratio

# Other libs
import random
import numpy as np
import galois

from functools import cached_property

class CMPR(FeedbackFunction):
    def __init__(self, config):
        #list of [0..sizes.. size]
        self.num_components = len(config)
        self.divisions = [] 

        shift_amount = 0
        self.fn_list = []

        for component in config[::-1]:
            # merge any CMPRs inside
            if isinstance(component,CMPR):
                self.divisions += [d + shift_amount for d in component.divisions[:-1]]
                self.num_components += component.num_components-1
                self.fn_list += [f.shift_indices(shift_amount) for f in component.fn_list]
            else: 
                self.divisions.append(shift_amount)
                self.fn_list += [XOR(f.shift_indices(shift_amount)) for f in component.fn_list]

            shift_amount += len(component.fn_list)
    
        self.size = len(self.fn_list)
        self.divisions.append(self.size)





    def generateChaining(self,
        template = None,
        modifier = None,
        **kwargs, 
    ):

        """
        parameters:

        block_probabilities
        input_densities
        input_minimums
        
        current_depth
        max_depth
        depth_mode

        component_distributions
        """

        # override the template used with any use entered variables.
        if template:
            combined_args = template(self)
            for k,v in kwargs.items():
                combined_args[k] = v
            kwargs = combined_args

        # pull these values out of the argument dictionary (not passed to random_function)
        input_densities = kwargs["input_densities"]
        input_minimums = kwargs["input_minimums"]
        del kwargs["input_densities"]
        del kwargs["input_minimums"]

        if "modifier" in kwargs:
            modifier = kwargs["modifier"]
            del kwargs["modifier"]

        # main loop
        for block_idx in range(1,self.num_components):
            has_chaining = []
            dist_used = min(len(combined_args["component_distributions"])-1, block_idx)
            block_component_dist = combined_args["component_distributions"][dist_used]


            # random pass
            for bit in self.blocks[block_idx]:
                # flip a coin to see if this bit is going to have input:
                if (random.random() <= input_densities[block_idx]):
                    has_chaining.append(bit)

                    # modify the args to add pass to random_function:
                    modified_args = {"allowed_blocks": self.blocks[:block_idx], **combined_args}
                    modified_args["component_distributions"] = block_component_dist
                    chaining_fn = random_function(**modified_args)

                    # allow some post generation deterministic modification
                    if modifier:
                        chaining_fn = modifier(chaining_fn, modified_args['allowed_blocks'])

                    self.fn_list[bit].add_arguments(chaining_fn)

            # minimum pass  
            if len(has_chaining) < input_minimums[block_idx]:
                inpt_bits = random.sample(
                    [bit for bit in self.blocks[block_idx] if not bit in has_chaining],
                    input_minimums[block_idx] - len(has_chaining)
                )

                for bit in inpt_bits:
                    # modify the args to add pass to random_function:
                    modified_args = {"allowed_blocks": self.blocks[:block_idx], **combined_args}
                    modified_args["component_distributions"] = block_component_dist

                    self.fn_list[bit].add_arguments(random_function(**modified_args))



    @property
    def has_chaining(self):
        return [len(f.args)-1 for f in self.fn_list]

    @property
    def component_feedback(self):
        return [f.args[0] for f in self.fn_list]

    @property
    def chaining_feedback(self):
        output = []
        for f in self.fn_list:
            if len(f.args) > 1:
                output.append(f.args[1])
            else:
                output.append(CONST(0))
        return output

    @cached_property 
    def blocks(self):
        block_list = []
        for d in range(len(self.divisions)-1):
            bits = list(range(self.divisions[d], self.divisions[d+1]))
            block_list.append(bits)
        return block_list[::-1]

    @property
    def update_matrices(self):
        matrices = []
        for b in range(len(self.blocks)):
            block = self.blocks[b]
            size = len(block)
            offset = self.divisions[-(b+1)]

            matrix = np.zeros([size,size], dtype = int)

            for inpt in block:
                for outpt in block:
                    #iterate through the VAR objects in the linear function portion
                    for leaf in self.component_feedback[outpt].inputs():
                        if leaf.index == inpt:
                            matrix[outpt-offset][inpt-offset] = 1
            
            matrices.append(matrix)
        return matrices

    @property
    def resolvent_matrices(self):
        resolvent_matrices = []

        unit = ResolventSolving.RationalPolynomial.unit()
        delay = ResolventSolving.RationalPolynomial(galois.Poly([1,0]),galois.Poly([1])) 
        for update_matrix in self.update_matrices:

            # Convert the update matrix to be over the Rational Polynomial Field
            converted_update_matrix = np.vectorize(ResolventSolving.RationalPolynomial.from_int)(update_matrix)
            converted_update_matrix.dtype = ResolventSolving.RationalPolynomial

            # Create the matrix to solve: (I xor UD) 
            field_matrix = np.asarray([delay]) * converted_update_matrix
            field_matrix += np.asarray([unit]) * ResolventSolving.field_eye(
                field = ResolventSolving.RationalPolynomial,
                size = update_matrix.shape[0]
            )
            
            # invert the matrix: (I xor UD)^{-1}
            resolvent_matrices.append(ResolventSolving.field_invert(
                field = ResolventSolving.RationalPolynomial,
                matrix = field_matrix
            ))

        return resolvent_matrices

    @property
    def propagation_matrices(self):
        propagation_matrices = []
        for resolvent_matrix in self.resolvent_matrices:
            mask = np.zeros_like(resolvent_matrix)
            mask[resolvent_matrix != 0] = 1

            propagation_matrices.append(mask)
        return propagation_matrices







    @cached_property
    def expected_period_ratio(self):
        numerator = 1
        denominator = 1
        already_seen = set()
        sizes = [self.divisions[i+1] - self.divisions[i]
                 for i in range(self.num_components)
                ]

        for s in sizes:
            denominator *= (2**(2*s))
            if s in already_seen: # this is wrong :/
                numerator *= (2**(s+1)-1)
            else:
                numerator *= (((2**s)-1)**2 + 1)
                already_seen.add(s)

        return numerator/denominator

    @cached_property
    def expected_period(self):
        numerator = 1
        denominator = 1
        already_seen = set()
        sizes = [self.divisions[i+1] - self.divisions[i]
                 for i in range(self.num_components)
                ]

        for s in sizes:
            denominator *= (2**s)
            if s in already_seen: # this is wrong :/
                numerator *= (2**(s+1)-1) 
            else:
                numerator *= (((2**s)-1)**2 + 1)
                already_seen.add(s)

        return numerator/denominator

    @cached_property
    def max_period(self):
        period = 1
        already_seen = set()
        sizes = [self.divisions[i+1] - self.divisions[i]
                 for i in range(self.num_components)
                ]
        
        for s in sizes:
            if s in already_seen:
                period *= 2
            else:
                already_seen.add(s)
                period *= (2**s-1) 
        
        return period










    def root_expressions(self, locked_list = None):
        expr_table = [RootExpression() for i in range(self.size)] # map: bit -> expression
        
        #fill in the following blocks:
        for block_idx in range(len(self.blocks)):
            block_RE = RootExpression()

            # map each bit to its corresponding block.
            for bit in self.blocks[block_idx]:
                
                modified_fn = self.fn_list[bit].remap_constants({
                    0: RootExpression.logical_zero(),
                    1: RootExpression.logical_one()
                })
                
                block_RE += modified_fn.eval_ANF(expr_table)

            #if this one isn't locked, extend it. 
            if (not locked_list) or (not locked_list[block_idx]):
                block_RE = block_RE.extend(JordanSet({len(self.blocks[block_idx]):1}, 1))
            
            #fill in table entries
            for bit in self.blocks[block_idx]:
                    expr_table[bit] = block_RE
        return expr_table




    def estimate_LC(self, output_bit, locked_list = None, benchmark = False):
        # locked-list is used to cancel effects of the locked registers.
        # the locked list contains the sizes of the locked MPRs

        from time import time_ns
        t1 = time_ns()
        REs = self.root_expressions(locked_list)
        bitRE = REs[output_bit]
        t2 = time_ns()

        if benchmark:
            print(f"Root Expression Generation: {len(bitRE.terms)} terms generated in {t2-t1} ns")

        t1 = time_ns()
        #get the length of the block bit is in.
        blocks = self.blocks
        for block in blocks:
            if output_bit in block:
                blockLen = len(block)
        
        #lower the minimum if this block is locked:
        if locked_list and blockLen in locked_list:
            blockLen = 1

        #add 1 to include the "True" / 1 in GF(2) value we ignore:
        upper = bitRE.upper(locked_list)

        #subtract safety_factor * size to account for some natural degeneracies:
        # TODO: need a better incorporation of the safety factor.

        lower = max(blockLen,bitRE.lower(locked_list))
        t2 = time_ns()
        if benchmark:
            print(f"Terms evaluated in {t2-t1} ns")

        return (lower,upper)




    @property
    def fixpoint(self):
        key = [0] * self.size
        for block_idx in range(self.num_components):
            # setup variables
            matrix_size = len(self.blocks[block_idx])
            update_matrix = self.update_matrices[block_idx]
            identity_matrix = np.eye(matrix_size, dtype = int)

            # Compute the chaining from known bits
            vector = [self.fn_list[i].eval(key) for i in self.blocks[block_idx]]

            # create matrices and invert
            galois_mat = galois.GF2(update_matrix) - galois.GF2(identity_matrix)
            galois_vect = galois.GF2(vector)
            sol = np.linalg.solve(galois_mat, galois_vect)

            # plug answer back into the key
            shift = self.blocks[block_idx][0]
            for bit in self.blocks[block_idx]:
                key[bit] = int(sol[bit - shift])
        return key



    def reverse_clock(self, state):
        key = [0] * self.size
        for block_idx in range(self.num_components):
            update_matrix = galois.GF2(self.update_matrices[block_idx])

            # Compute the chaining from known bits
            chaining_vector = galois.GF2([self.fn_list[i].eval(key) for i in self.blocks[block_idx]])
            target_vector = galois.GF2([state[i] for i in self.blocks[block_idx]])

            # solve the system
            sol = np.linalg.solve(update_matrix, target_vector - chaining_vector)

            # plug answer back into the key
            shift = self.blocks[block_idx][0]
            for bit in self.blocks[block_idx]:
                key[bit] = int(sol[bit - shift])
        return key





    #writes a VHDL file (special formatting for CMPRs)
    #Credit: Anna Hemingway
    def write_VHDL(self, filename):
        with open(filename, "w") as f:
            f.write(f"""
library ieee;
use ieee.std_logic_1164.all;

entity fpr is
    port (
    i_clk :in std_logic;
    i_rst : in std_logic;
    i_seed_data: in std_logic_vector( {self.size - 1} downto 0);
    output: out std_logic_vector({self.size - 1} downto 0)
    );
end entity fpr;

architecture run of fpr is

    signal currstate, nextstate:std_logic_vector({self.size - 1} downto 0);


begin

    statereg: process(i_clk, i_rst)
    begin
        if (i_rst = '1') then
            currstate <= i_seed_data;
        elsif (i_clk = '1' and i_clk'event) then
            currstate <= nextstate;
        end if;
    end process;\n""")
        
            for i in range(self.size - 1, -1 , -1):

                if not self.has_chaining[i]:
                    f.write(f"    nextstate({str(i)}) <= {self.component_feedback[i].generate_VHDL()};\n")
                else:
                    f.write(f"    nextstate({str(i)}) <= {self.component_feedback[i].generate_VHDL()} XOR \n" +
                            f"        {self.chaining_feedback[i].generate_VHDL()};\n")
            f.write("""
    output <= currstate;

end run;

""")