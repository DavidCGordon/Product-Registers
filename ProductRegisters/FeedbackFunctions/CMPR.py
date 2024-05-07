# Documentation:
from typing import Iterable

# Simulation Boilerplate
from ProductRegisters.BooleanLogic import ANF_spec_repr, XOR, AND, CONST, VAR
from ProductRegisters.FeedbackFunctions import FeedbackFunction
from ProductRegisters.FeedbackFunctions import MPR

# Linear Complexity and Monomial estimation
from ProductRegisters.Tools.RootCounting.MonomialProfile import TermSet,MonomialProfile
from ProductRegisters.Tools.RootCounting.JordanSet import JordanSet
from ProductRegisters.Tools.RootCounting.RootExpression import RootExpression

# Other analysis
import ProductRegisters.Tools.ResolventSolving as ResolventSolving
from ProductRegisters.Tools.MersenneTools import expected_period, expected_period_ratio

# Other libs
import random
import numpy as np
import galois as gl
import time

from functools import cached_property

class CMPR(FeedbackFunction):
    def __init__(self, 
        components: Iterable[MPR]
    ):
        """
        Constructor for the CMPR family of Feedbacl Functions

        Args:
          components: 
            An iterable of MPR objects to be used as the components of the CMPR. No error will be raised
            if feedback functions which are not MPRs are used, and this can sometimes be a useful hack.
            however not all methods are guaranteed to work if non-MPR components are used
        Returns:
            None
        """
        self.num_components = len(components)
        self.divisions = [] 

        shift_amount = 0
        self.fn_list = []

        for component in components[::-1]:
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





    def generateChaining(self,template):
        chaining_logic = template(self)

        for bit, fn in chaining_logic.items():
            self.fn_list[bit].add_arguments(fn)


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
                output.append(XOR(*f.args[1:]))
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

    @cached_property
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

    @cached_property
    def resolvent_matrices(self):
        resolvent_matrices = []
        for update_matrix in self.update_matrices:

            # Convert the update matrix to be over the Rational Polynomial Field
            converted_update_matrix = np.vectorize(ResolventSolving.SequenceTransform.from_int)(update_matrix)
            converted_update_matrix.dtype = ResolventSolving.SequenceTransform

            # (I xor UD)^{-1}):
            # meant to be multiplied by (DC(D) xor B[0])
            unit = ResolventSolving.SequenceTransform.one()
            delay = ResolventSolving.SequenceTransform.delay()
            
            field_matrix = np.asarray([delay]) * converted_update_matrix
            field_matrix += np.asarray([unit]) * ResolventSolving.field_eye(
                field = ResolventSolving.SequenceTransform,
                size = update_matrix.shape[0],
            )

            # invert the matrix and append
            resolvent_matrices.append(ResolventSolving.field_invert(
                field = ResolventSolving.SequenceTransform,
                matrix = field_matrix
            ))

        return resolvent_matrices

    @cached_property
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










    def monomial_profiles(self):
        prof_table = [MonomialProfile() for i in range(self.size)] # map: bit -> expression
        block_table = [MonomialProfile() for i in range(len(self.blocks))]
        
        #fill in the following blocks:
        for block_id in range(len(self.blocks)):
            monomial_profile = MonomialProfile.from_merged(
                fn_list = [self.fn_list[i] for i in self.blocks[block_id]],
                blocks = self.blocks
            )

            block_fn = monomial_profile.to_BooleanFunction().remap_constants({
                0: MonomialProfile.logical_zero(),
                1: MonomialProfile.logical_one()
            })
                
            block_table[block_id] = block_fn.eval_ANF(block_table)

            # if this one isn't locked, extend it. 
            block_table[block_id] += MonomialProfile([TermSet(
                {block_id: len(self.blocks[block_id])},
                {block_id: 1}
            )])
            
            #fill in table entries
            for bit in self.blocks[block_id]:
                prof_table[bit] = block_table[block_id].__copy__()
        return prof_table


    def root_expressions(self, locked_list = None, verbose = False):
        expr_table = [RootExpression() for i in range(self.size)] # map: bit -> expression
        block_table = [RootExpression() for i in range(len(self.blocks))]
        
        #fill in the following blocks:
        for block_id in range(len(self.blocks)):
            start_time = time.time()
            
            if verbose: print("Starting Monomial Profiling")

            monomial_profile = MonomialProfile.from_merged(
                fn_list = [self.fn_list[i] for i in self.blocks[block_id]],
                blocks = self.blocks
            )

            block_fn = monomial_profile.to_BooleanFunction().remap_constants({
                0: RootExpression.logical_zero(),
                1: RootExpression.logical_one()
            })

            if verbose:
                print(f"Monomial Profiled: {block_fn.dense_str()}")
                print(f"Monomial Profiling Time: {time.time()-start_time}\n")

            start_time = time.time()

            if verbose:
                print("Starting ANF Composition")    

            block_table[block_id] = block_fn.eval_ANF(block_table)

            # if this one isn't locked, extend it. 
            if (not locked_list) or (not locked_list[block_id]):
                block_table[block_id] = block_table[block_id].extend(
                    JordanSet({len(self.blocks[block_id]):1}, 1)
                )

            #fill in table entries
            for bit in self.blocks[block_id]:
                expr_table[bit] = block_table[block_id].__copy__()
            
            if verbose:
                print(f'Block {block_id} finished  -  Num Terms: {len(block_table[block_id].terms)}')
                print(f"ANF Composition Time: {time.time()-start_time}\n\n\n")

        return expr_table
    

    def estimate_LC(self, output_bit, locked_list = None, verbose = False):
        # locked-list is used to cancel effects of the locked registers.
        # the locked list contains the sizes of the locked MPRs

        from time import time_ns
        t1 = time_ns()
        REs = self.root_expressions(locked_list,verbose = verbose)
        bitRE = REs[output_bit]
        t2 = time_ns()

        if verbose:
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

        upper = bitRE.upper(locked_list)
        lower = max(blockLen,bitRE.lower(locked_list))

        t2 = time_ns()
        if verbose:
            print(f"Terms evaluated in {t2-t1} ns")

        return (lower,upper)









    @property
    def fixpoint(self):
        fixed_state = [0] * self.size
        for block_idx in range(self.num_components):
            # compute the matrix (U-I)
            matrix_size = len(self.blocks[block_idx])
            update_matrix = gl.GF2(self.update_matrices[block_idx])
            identity_matrix = gl.GF2(np.eye(matrix_size, dtype = int))
            difference_matrix = update_matrix - identity_matrix

            # compute the target chaining vector from known bits
            chaining_vector = gl.GF2([self.fn_list[i].eval(fixed_state) for i in self.blocks[block_idx]])

            # solve the system (U-I)x = C.
            # expanding gives: Ux + x = C.
            # swapping terms:  Ux + C = x
            sol = np.linalg.solve(difference_matrix, chaining_vector)

            # write answer back into the fixed state vector
            shift = self.blocks[block_idx][0]
            for bit in self.blocks[block_idx]:
                fixed_state[bit] = int(sol[bit - shift])
        return fixed_state


    def reverse_clock(self, state):
        prev_state = [0] * self.size
        for block_idx in range(self.num_components):
            update_matrix = gl.GF2(self.update_matrices[block_idx])

            # Compute the chaining from known bits
            chaining_vector = gl.GF2([self.fn_list[i].eval(prev_state) for i in self.blocks[block_idx]])
            target_vector = gl.GF2([state[i] for i in self.blocks[block_idx]])

            # solve the system Ux = T - C
            # rearranging:     Ux + C = T 
            sol = np.linalg.solve(update_matrix, target_vector - chaining_vector)

            # plug answer back into the key
            shift = self.blocks[block_idx][0]
            for bit in self.blocks[block_idx]:
                prev_state[bit] = int(sol[bit - shift])
        return prev_state





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