from __future__ import annotations
from typing import Any, Iterable

from PyPR.BooleanLogic import BooleanANF, XOR, AND, CONST, VAR
from PyPR.FeedbackFunctions import FeedbackFunction
from PyPR.FeedbackFunctions import MPR

from PyPR.Tools.RootCounting.MonomialProfile import TermSet, MonomialProfile
from PyPR.Tools.RootCounting.JordanSet import JordanSet
from PyPR.Tools.RootCounting.RootExpression import RootExpression
import PyPR.Tools.RootCounting.MeshOptimization as mesh_optimization

import PyPR.Tools.ResolventSolving as ResolventSolving
from PyPR.BooleanLogic.BooleanGF import BooleanGF
from PyPR.Tools.MersenneTools import expected_period, expected_period_ratio, max_period, cycle_lengths

import random
import numpy as np
import galois as gl
import time

from functools import cached_property

class CMPR(FeedbackFunction):
    """A Composite Mersenne Product Register.

    A CMPR couples several MPR components via nonlinear chaining logic.
    Each component operates on its own block of bits as an independent
    field multiplier (see :class:`MPR`); chaining functions inject
    nonlinear dependencies from higher-index blocks into lower-index
    blocks, breaking the linearity of any single component and increasing
    the linear complexity of the output sequence.

    The blocks are ordered so that block 0 is the highest-index
    (rightmost) block and chaining propagates from higher-index blocks
    toward lower ones.

    :ivar size: Total number of bits across all component blocks.
    :vartype size: int
    :ivar num_components: The number of MPR component blocks.
    :vartype num_components: int
    :ivar primitive_polynomials: The primitive polynomial of each component
        block, indexed by block number.
    :vartype primitive_polynomials: list[list[int] | None]
    :ivar update_polynomials: The update polynomial of each component
        block, indexed by block number.
    :vartype update_polynomials: list[list[int] | None]
    :ivar divisions: Internal list of bit-index boundaries between blocks.
    :vartype divisions: list[int]
    """

    num_components: int
    primitive_polynomials: list[list[int] | None]
    update_polynomials: list[list[int] | None]
    divisions: list[int]

    def __init__(self,
        components: list["FeedbackFunction | MPR | CMPR"]
    ) -> None:
        """Construct a CMPR from a list of MPR (or nested CMPR) components.

        Components are provided in high-to-low order (the first element
        becomes the highest-index block). Nested CMPRs are flattened
        automatically. A non-MPR feedback function may be supplied as the
        first (highest-index) component only this is occasionally useful as a
        hack, but not all methods will work with non-MPR components.

        :param components: The MPR (or CMPR) components to compose.
        :type components: list[MPR | CMPR]
        """
        self.num_components = len(components)

        self.primitive_polynomials = []
        self.update_polynomials = []
        self.divisions = [] 

        shift_amount = 0
        self.fn_list = []
        for i,component in enumerate(components[::-1]):
            # merge any CMPRs inside
            if isinstance(component,CMPR):
                self.divisions += [d + shift_amount for d in component.divisions[:-1]]
                self.primitive_polynomials += component.primitive_polynomials
                self.update_polynomials += component.update_polynomials
                self.num_components += component.num_components-1
                self.fn_list += [f.shift_indices(shift_amount) for f in component.fn_list]
            elif isinstance(component,MPR): 
                self.divisions.append(shift_amount)
                self.primitive_polynomials += [component.primitive_polynomial]
                self.update_polynomials += [component.update_polynomial]
                self.fn_list += [XOR(f.shift_indices(shift_amount)) for f in component.fn_list]
            else:
                # only allow different types for the first component
                if i == len(components)-1:
                    self.divisions.append(shift_amount)
                    self.primitive_polynomials += [None]
                    self.update_polynomials += [None]
                    self.fn_list += [XOR(f.shift_indices(shift_amount)) for f in component.fn_list]
                else:
                    raise TypeError(
                        f"Other than the component at index 0, all CMPR components must be "
                        f"either MPRs or CMPRs, not {type(component)}"
                    )
            shift_amount += len(component.fn_list)

        self.size = len(self.fn_list)
        self.divisions.append(self.size)
        
        # reverse polynomials to match block indices
        self.primitive_polynomials = self.primitive_polynomials[::-1]
        self.update_polynomials = self.update_polynomials[::-1]

    def generateChaining(self,
        template: Any
    ) -> None:
        """Install chaining logic from a template function.

        The template is a callable that receives this CMPR and returns a
        dict mapping bit indices to boolean functions. Each returned
        function is added as an additional argument to the corresponding
        bit's existing feedback, introducing cross-block dependencies.

        :param template: A chaining template callable. See
            :mod:`PyPR.BooleanLogic.ChainingGeneration` for built-in
            templates.
        :type template: Callable[[CMPR], dict[int, BooleanFunction]]
        """
        chaining_logic = template(self)

        for bit, fn in chaining_logic.items():
            self.fn_list[bit].add_arguments(fn)

    def update_MPR(self,
        mpr_index: int,
        new_update_poly: list[int]
    ) -> None:
        """Replace the update polynomial of a single component block.

        Rebuilds the component's feedback functions from the new update
        polynomial while preserving any chaining logic already installed.
        Cached matrices (update, resolvent, propagation) are invalidated.

        :param mpr_index: The block index of the component to update.
        :type mpr_index: int
        :param new_update_poly: The new update polynomial, of length equal
            to the block size.
        :type new_update_poly: list[int]
        :raises ValueError: If the block at `mpr_index` is not an MPR, or
            if `new_update_poly` has the wrong length.
        """
        if mpr_index == 0 and self.primitive_polynomials[0] == None:
            raise ValueError(f"Can't update Update Polynomial for block 0, because it is not an MPR.")
        if len(new_update_poly) != len(self.blocks[mpr_index]):
            raise ValueError(f"Update polynomial must be {len(self.blocks[mpr_index])} bits.")
        
        # create new MPR
        new_mpr = MPR(
            len(self.blocks[mpr_index]),
            self.primitive_polynomials[mpr_index], #type: ignore (None case handled above)
            new_update_poly
        )

        # update CMPR functions
        shift = 2
        cmpr_bits = self.blocks[mpr_index]
        for cmpr_bit, mpr_bit in zip(cmpr_bits, range(len(cmpr_bits))):
            self.fn_list[cmpr_bit] = XOR(
                new_mpr.fn_list[mpr_bit].shift_indices(shift),
                *self.fn_list[cmpr_bit].args[1:]
            )
        
        # update stored polynomials
        self.update_polynomials[mpr_index] = new_update_poly

        # refresh cached properties
        if 'update_matrices' in self.__dict__: del self.update_matrices
        if 'resolvent_matrices' in self.__dict__: del self.resolvent_matrices
        if 'propagation_matrices' in self.__dict__: del self.propagation_matrices

    @property
    def has_chaining(self) -> list[int]:
        """The number of chaining terms on each bit (0 if unchained).

        :return: A per-bit list of chaining-term counts.
        :rtype: list[int]
        """
        return [len(f.args)-1 for f in self.fn_list]

    @property
    def component_feedback(self) -> list[Any]:
        """The linear (MPR) portion of each bit's feedback function.

        :return: A per-bit list of the component-only feedback (the first
            argument of each bit's top-level XOR).
        :rtype: list[BooleanFunction]
        """
        return [f.args[0] for f in self.fn_list]

    @property
    def chaining_feedback(self) -> list[Any]:
        """The nonlinear chaining portion of each bit's feedback function.

        For bits with no chaining, returns ``CONST(0)``.

        :return: A per-bit list of the chaining-only feedback.
        :rtype: list[BooleanFunction]
        """
        output = []
        for f in self.fn_list:
            if len(f.args) > 1:
                output.append(XOR(*f.args[1:]))
            else:
                output.append(CONST(0))
        return output

    @cached_property
    def blocks(self) -> list[list[int]]:
        """The bit-index ranges of each component block.

        Block 0 is the highest-index (rightmost) block. Each entry is a
        list of the bit indices belonging to that block.

        :return: A list of block-index lists, one per component.
        :rtype: list[list[int]]
        """
        block_list = []
        for d in range(len(self.divisions)-1):
            bits = list(range(self.divisions[d], self.divisions[d+1]))
            block_list.append(bits)
        return block_list[::-1]

    @cached_property
    def update_matrices(self) -> list[np.ndarray[tuple[int],np.dtype[np.uint8]]]:
        """The GF(2) update matrix of each component block.

        For block b, the matrix U satisfies ``next_block = U @ block``
        (over GF(2)) on the linear (component-only) portion of the
        feedback. Entry ``U[i, j]`` is 1 iff bit i's component feedback
        references bit j as a VAR leaf.

        :return: A list of square GF(2) matrices, one per block.
        :rtype: list[numpy.ndarray]
        """
        matrices = []
        for b in range(len(self.blocks)):
            block = self.blocks[b]
            size = len(block)
            offset = self.divisions[-(b+1)]

            matrix = np.zeros([size,size], dtype = np.uint8)

            for inpt in block:
                for outpt in block:
                    #iterate through the VAR objects in the linear function portion
                    for leaf in self.component_feedback[outpt].inputs():
                        if leaf.index == inpt:
                            matrix[outpt-offset][inpt-offset] = 1

            matrices.append(matrix)
        return matrices

    @cached_property
    def resolvent_matrices(self) -> list[np.ndarray[tuple[int],np.dtype[np.object_]]]:
        """The resolvent matrix (I + UD)^{-1} for each component block.

        Computed over the formal power series ring GF(2)[[D]], where D is
        the delay operator. The resolvent captures how an input
        perturbation at a given bit propagates forward through time within
        its block.

        :return: A list of resolvent matrices over GF(2)[[D]], one per
            block.
        :rtype: list[numpy.ndarray]
        """
        resolvent_matrices = []
        for update_matrix in self.update_matrices:

            # Convert the update matrix to be over the Rational Polynomial Field
            converted_update_matrix = np.vectorize(ResolventSolving.BooleanGF.from_int)(update_matrix)
            converted_update_matrix.dtype = ResolventSolving.BooleanGF

            # (I xor UD)^{-1}):
            # meant to be multiplied by (DC(D) xor B[0])
            unit = ResolventSolving.BooleanGF.one()
            delay = ResolventSolving.BooleanGF.delay()

            field_matrix = np.asarray([delay]) * converted_update_matrix
            field_matrix += np.asarray([unit]) * ResolventSolving.field_eye(
                field = ResolventSolving.BooleanGF,
                size = update_matrix.shape[0],
            )

            # invert the matrix and append
            resolvent_matrices.append(ResolventSolving.field_invert(
                field = ResolventSolving.BooleanGF,
                matrix = field_matrix
            ))

        return resolvent_matrices

    @cached_property
    def propagation_matrices(self) -> list[np.ndarray]:
        """The support mask of each resolvent matrix.

        Entry (i, j) is 1 if the corresponding resolvent entry is
        nonzero, indicating that bit j's perturbation eventually affects
        bit i within the same block.

        :return: A list of binary masks, one per block.
        :rtype: list[numpy.ndarray]
        """
        propagation_matrices = []
        for resolvent_matrix in self.resolvent_matrices:
            mask = np.zeros_like(resolvent_matrix) # ignoring error cause by dtype being object 
            mask[resolvent_matrix != resolvent_matrix.dtype.zero] = 1 #type: ignore
            propagation_matrices.append(mask)
        return propagation_matrices




    @cached_property
    def expected_period_ratio(self) -> float:
        """The expected fraction of the maximum period achieved by a random
        chaining configuration with these block sizes.

        :return: A probability in [0, 1].
        :rtype: float
        """
        sizes = [len(block) for block in self.blocks]
        return expected_period_ratio(sizes)

    @cached_property
    def expected_period(self) -> float:
        """The expected period for a random chaining configuration with
        these block sizes.

        :return: The expected period.
        :rtype: int
        """
        sizes = [len(block) for block in self.blocks]
        return expected_period(sizes)

    @cached_property
    def max_period(self) -> int:
        """The maximum achievable period for these block sizes.

        Equal to the LCM of (2^{n_i} - 1) across all component blocks.

        :return: The maximum period.
        :rtype: int
        """
        sizes = [len(block) for block in self.blocks]
        return max_period(sizes)

    @cached_property
    def cycle_lengths(self) -> list[tuple[int, int]]:
        """The possible cycle lengths and their multiplicities for these
        block sizes.

        :return: A list of (cycle_length, count) pairs.
        :rtype: list[tuple[int, int]]
        """
        sizes = [len(block) for block in self.blocks]
        return cycle_lengths(sizes)






    def monomial_profiles(self,
        verbose: bool = False,
        force_default: bool = False
    ) -> list[MonomialProfile]:
        """Compute a monomial profile for each bit of the register.

        The monomial profile tracks which monomials (products of roots
        from distinct component blocks) can appear in the exponential
        representation of each bit's output sequence. This provides an
        upper bound on the algebraic structure and, in turn, the linear
        complexity.

        When the block configuration is amenable (no 1-bit blocks, no
        repeated sizes, simple chaining), a mesh optimization is used for
        faster computation. Otherwise falls back to the default
        composition-based algorithm.

        :param verbose: If True, print progress information.
        :type verbose: bool
        :param force_default: If True, skip the mesh optimization even
            when it would be valid.
        :type force_default: bool
        :return: A per-bit list of monomial profiles.
        :rtype: list[MonomialProfile]
        """
        block_sizes = set()
        use_mesh_optimization = True

        for block_id in range(len(self.blocks)):
            # check if optimization not valid due to 1 bit MPR:
            if len(self.blocks[block_id]) == 1:
                use_mesh_optimization = False
                if verbose: print("Found 1-bit MPR")
                break

            # check if optimization not valid due to repeated sizes
            if len(self.blocks[block_id]) in block_sizes:
                use_mesh_optimization = False
                if verbose: print("Found duplicate size")
                break
            block_sizes.add(len(self.blocks[block_id]))
            
            # check if optimization not valid due to complex chaining
            allowed_bits = set(self.blocks[block_id]) | set(self.blocks[block_id-1])
            used_bits = set().union(*(
                self.fn_list[bit].idxs_used() for bit in self.blocks[block_id]
            ))
            
            if not (used_bits <= allowed_bits):
                use_mesh_optimization = False
                if verbose: print("Found non-simple chaining function")
                break

        if use_mesh_optimization and not force_default:
            return self._mp_mesh_optimization(verbose)
        else:
            return self._mp_default(verbose)

    def _mp_default(self, verbose: bool = False) -> list[MonomialProfile]:
        if verbose: print("Running default monomial profile algorithm")
        prof_table = [MonomialProfile() for i in range(self.size)] # map: bit -> expression
        block_table = [MonomialProfile() for i in range(len(self.blocks))]
        
        #fill in the following blocks:
        for block_id in range(len(self.blocks)):
            start_time = time.time()
            if verbose: print("Profiling Chaining")

            chaining_profile = MonomialProfile.from_merged(
                fn_list = [self.fn_list[i] for i in self.blocks[block_id]],
                blocks = self.blocks
            ).to_BooleanFunction()

            if verbose:
                print(f"Chaining Profile: {chaining_profile.dense_str()}")
                print(f"Profiling Time: {time.time()-start_time}\n")

            # combine function
            block_fn = chaining_profile.remap_constants([
                (0, MonomialProfile.logical_zero()),
                (1, MonomialProfile.logical_one())
            ])

            start_time = time.time()
            if verbose:
                print("Starting ANF Composition")    
            
            # compose MPs into the block table:
            block_table[block_id] = block_fn.eval_ANF(block_table) 
            block_table[block_id] += MonomialProfile([TermSet(
                {block_id: len(self.blocks[block_id])},
                {block_id: 1}
            )])
            
            # fill in table entries
            for bit in self.blocks[block_id]:
                prof_table[bit] = block_table[block_id].__copy__()

            if verbose:
                num_terms = len(block_table[block_id].terms)
                print(f'Block {block_id} finished  -  Num Terms: {num_terms}')
                print(f"ANF Composition Time: {time.time()-start_time}\n\n\n")

        return prof_table

    def _mp_mesh_optimization(self, verbose: bool = False) -> list[Any]:
        if verbose: print("Running monomial profile algorithm with the mesh optimization")
                
        expr_table: list[Any] = [None for i in range(self.size)]

        # compute degree list used for the pass:
        degrees = []
        sizes = [len(block) for block in self.blocks]
        constants_possible = [False for block in self.blocks]
        for block_id in range(len(self.blocks)):
            chaining_profile = MonomialProfile.from_merged(
                fn_list = [self.fn_list[i] for i in self.blocks[block_id]],
                blocks = self.blocks
            ).to_BooleanFunction()

            # Repeated Vars => can't use normal degree
            degree = 0
            for term in chaining_profile.args:
                if type(term) == VAR:
                    degree = max(degree, 1)
                elif type(term) == CONST and term.value == 1:
                    constants_possible[block_id] = True
                else:
                    degree = max(degree,len(term.args))
            degrees.append(degree)
        

        # write the first block directly:
        size = len(self.blocks[0])
        for bit in self.blocks[0]:
            expr_table[bit] = MonomialProfile([TermSet({0:size},{0:1})])

        # calculate the RE for each subsequent block:
        for block_id in range(1,len(self.blocks)):
            start_time = time.time()

            if verbose:
                print(f"Iterating over mesh for block {block_id} (degree {degrees[block_id]})")

            monomial_profile = mesh_optimization.mp_compute_single_mesh(
                sizes[:block_id+1],
                degrees[:block_id+1],
            )

            if constants_possible[block_id]:
                monomial_profile += MonomialProfile.logical_one()

            end_time = time.time()

            # write to table
            for bit in self.blocks[block_id]:
                expr_table[bit] = monomial_profile.__copy__()

            if verbose:
                num_terms = len(monomial_profile.terms)
                print(f'Block {block_id} finished  -  Num Terms: {num_terms}')
                print(f"ANF Composition Time: {end_time-start_time}")
                print(f"Copying Time: {time.time()-end_time}\n\n\n")

        return expr_table









    def root_expressions(self,
        locked_list: list[int] | None = None,
        verbose: bool = False,
        force_default: bool = False
    ) -> list[RootExpression]:
        """Compute a root expression for each bit of the register.

        The root expression bounds which roots can appear in the
        exponential (Binet-style) representation of each bit's output
        sequence over GF(2). This provides upper and lower bounds on the
        linear complexity, which measures the shortest LFSR that
        reproduces the sequence.

        Blocks listed in `locked_list` are treated as fixed (their root
        contributions are not extended), useful for analyzing the effect of
        chaining with some components held constant.

        :param locked_list: A list of flags (one per block), where a
            truthy value means the block is unlocked (contributes roots).
            If None, all blocks are unlocked.
        :type locked_list: list[int] | None
        :param verbose: If True, print progress information.
        :type verbose: bool
        :param force_default: If True, skip the mesh optimization even
            when it would be valid.
        :type force_default: bool
        :return: A per-bit list of root expressions.
        :rtype: list[RootExpression]
        """
        block_sizes = set()
        use_mesh_optimization = True

        for block_id in range(len(self.blocks)):
            # check if optimization not valid due to 1 bit MPR:
            if len(self.blocks[block_id]) == 1:
                use_mesh_optimization = False
                if verbose: print("Found 1-bit MPR")
                break

            # check if optimization not valid due to repeated sizes
            if len(self.blocks[block_id]) in block_sizes:
                use_mesh_optimization = False
                if verbose: print("Found duplicate size")
                break
            block_sizes.add(len(self.blocks[block_id]))
            
            # check if optimization not valid due to complex chaining
            allowed_bits = set(self.blocks[block_id]) | set(self.blocks[block_id-1])
            used_bits = set().union(*(
                self.fn_list[bit].idxs_used() for bit in self.blocks[block_id]
            ))
            
            if not (used_bits <= allowed_bits):
                use_mesh_optimization = False
                if verbose: print("Found non-simple chaining function")
                break

        if use_mesh_optimization and not force_default:
            return self._re_mesh_optimization(locked_list, verbose)
        else:
            return self._re_default(locked_list, verbose)

    def _re_default(self,
        locked_list: list[int] | None = None,
        verbose: bool = False
    ) -> list[RootExpression]:
        if verbose: print("Running default root expression algorithm")
        expr_table = [RootExpression({}) for i in range(self.size)] # map: bit -> expression
        block_table = [RootExpression({}) for i in range(len(self.blocks))]
        
        #fill in the following blocks:
        for block_id in range(len(self.blocks)):
            start_time = time.time()
            
            if verbose: print("Profiling Chaining")

            monomial_profile = MonomialProfile.from_merged(
                fn_list = [self.fn_list[i] for i in self.blocks[block_id]],
                blocks = self.blocks
            ).to_BooleanFunction()

            if verbose:
                print(f"Chaining Profile: {monomial_profile.dense_str()}")
                print(f"Profiling Time: {time.time()-start_time}\n")

            block_fn = monomial_profile.remap_constants([
                (0, RootExpression.logical_zero()),
                (1, RootExpression.logical_one())
            ])

            start_time = time.time()

            if verbose:
                print("Starting ANF Composition")    

            block_table[block_id] = block_fn.eval_ANF(block_table)

            # if this one isn't locked, extend it. 
            if (not locked_list) or (locked_list[block_id]):
                size = len(self.blocks[block_id])

                # size 1 blocks are the same as constant/logical 1
                # and we use empty notation for true to avoid clutter
                if size == 1:
                    block_table[block_id] = block_table[block_id].extend(
                        JordanSet({},{1})
                    )
                else:
                    block_table[block_id] = block_table[block_id].extend(
                        JordanSet({size:1},[1])
                    )

            #fill in table entries
            for bit in self.blocks[block_id]:
                expr_table[bit] = block_table[block_id].__copy__()
            
            if verbose:
                num_terms = sum(len(table_entry) for table_entry in block_table[block_id].root_table.values())
                print(f'Block {block_id} finished  -  Num Terms: {num_terms}')
                print(f"ANF Composition Time: {time.time()-start_time}\n\n\n")

        return expr_table
    
    def _re_mesh_optimization(self,
        locked_list: list[int] | None = None,
        verbose: bool = False
    ) -> list[Any]:
        if verbose: print("Running root expression algorithm with the mesh optimization")
                
        expr_table: list[Any] = [None for i in range(self.size)]

        # compute lists used for the pass:
        degrees = []
        sizes = [len(block) for block in self.blocks]
        constants_possible = [False for block in self.blocks]
        for block_id in range(len(self.blocks)):
            chaining_profile = MonomialProfile.from_merged(
                fn_list = [self.fn_list[i] for i in self.blocks[block_id]],
                blocks = self.blocks
            ).to_BooleanFunction()

            block_fn = chaining_profile
            degree = 0
            for term in block_fn.args:
                if hasattr(term,'args'):
                    degree = max(degree,len(term.args))
                elif type(term) == CONST and term.value == 1:
                    constants_possible[block_id] = True
            degrees.append(degree)
        
        # if no locked list is passed, default to all registers unlocked:
        if not locked_list:
            locked_list = [1]*len(sizes)

        # write the first block directly:
        size = len(self.blocks[0])
        for bit in self.blocks[0]:
            js = JordanSet({size:1}, set([1]))
            expr_table[bit] = RootExpression({(size,):set([js])})

        # calculate the RE for each subsequent block:
        for block_id in range(1,len(self.blocks)):
            start_time = time.time()

            if verbose:
                print(f"Iterating over mesh for block {block_id} (degree {degrees[block_id]})")

            root_expression = mesh_optimization.re_compute_single_mesh(
                sizes[:block_id+1],
                degrees[:block_id+1],
                locked_list[:block_id+1]
            )

            if constants_possible[block_id]:
                root_expression += RootExpression.logical_one()

            end_time = time.time()

            # write to table
            for bit in self.blocks[block_id]:
                expr_table[bit] = root_expression.__copy__()

            if verbose:
                num_terms = sum(len(table_entry) for table_entry in root_expression.root_table.values())
                print(f'Block {block_id} finished  -  Num Terms: {num_terms}')
                print(f"ANF Composition Time: {end_time-start_time}")
                print(f"Copying Time: {time.time()-end_time}\n\n\n")

        return expr_table

               

    def estimate_LC(self,
        output_bit: int,
        locked_list: list[int] | None = None,
        verbose: bool = False
    ) -> tuple[int, int]:
        """Estimate the linear complexity bounds for a single output bit.

        Computes root expressions for the register and evaluates lower and
        upper bounds on the linear complexity of the specified bit's
        output sequence. Locked blocks (via `locked_list`) have their root
        contributions suppressed, which tightens the bound to reflect only
        the unlocked components' contributions.

        :param output_bit: The bit index to analyze.
        :type output_bit: int
        :param locked_list: A list of flags (one per block) indicating
            which blocks are unlocked. If None, all blocks are unlocked.
        :type locked_list: list[int] | None
        :param verbose: If True, print timing information.
        :type verbose: bool
        :return: A (lower, upper) pair bounding the linear complexity.
        :rtype: tuple[int, int]
        """
        t1 = time.time()
        REs = self.root_expressions(locked_list,verbose = verbose)
        bitRE = REs[output_bit]
        t2 = time.time()

        if verbose:
            print(f"Total RE generation time: {t2-t1} s\n\n\n")

        t1 = time.time()
        #get the length of the block bit is in.
        blocks = self.blocks
        for block in blocks:
            if output_bit in block:
                blockLen = len(block)
        
        #lower the minimum if this block is locked:
        if locked_list and blockLen in locked_list:
            blockLen = 1
                  
        upper = bitRE.upper()
        lower = max(blockLen,bitRE.lower())

        t2 = time.time()
        if verbose:
            print(f"Terms evaluated in {t2-t1} s")

        return (lower,upper)









    @property
    def fixpoint(self) -> list[int]:
        """The unique fixed point of the register (the state s such that
        F(s) = s).

        Solves (U - I)x = C block by block, where U is the update matrix
        and C is the chaining contribution from already-solved blocks.
        Requires that (U - I) is invertible for each block, which holds
        when the update polynomial has no fixed points.

        :return: The fixed-point state vector.
        :rtype: list[int]
        """
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


    def reverse_clock(self,
        state: list[int]
    ) -> list[int]:
        """Compute the predecessor state (the state that clocks into the
        given one).

        Solves Ux = target - chaining block by block, where U is the
        update matrix and the chaining contribution is computed from
        already-solved blocks.

        :param state: The current state vector.
        :type state: list[int]
        :return: The predecessor state vector.
        :rtype: list[int]
        """
        prev_state = [0] * self.size
        for block_idx in range(self.num_components):
            update_matrix = gl.GF2(self.update_matrices[block_idx])

            # Compute the chaining from known bits
            chaining_vector = gl.GF2([self.fn_list[i].eval(prev_state) for i in self.blocks[block_idx]])
            target_vector = gl.GF2([state[i] for i in self.blocks[block_idx]])

            # solve the system Ux = Target - Chaining
            # rearranging:     Ux + Chaining = Target
            sol = np.linalg.solve(update_matrix, target_vector - chaining_vector)

            # plug answer back into the key
            shift = self.blocks[block_idx][0]
            for bit in self.blocks[block_idx]:
                prev_state[bit] = int(sol[bit - shift])
        return prev_state


    def write_VHDL(self,
        filename: str,
        include_mpr: bool = True
    ) -> None:
        """Write a VHDL entity implementing this CMPR.

        Generates a synthesizable VHDL file with a clocked state register
        and combinational next-state logic derived from the component and
        chaining feedback functions. Credit: Anna Hemingway.

        :param filename: The output file path.
        :type filename: str
        :param include_mpr: If True, include the MPR (component) feedback
            in the output. If False, emit only the chaining logic.
        :type include_mpr: bool
        """
        overrides = {}
        vhdl_str = "\n    "
        for i in range(self.size - 1, -1 , -1):
            
            # write the current function:
            if self.has_chaining[i] and include_mpr:
                chaining_lines = self.chaining_feedback[i].merge_redundant().generate_VHDL(
                    output_name = f"next_state({i})",
                    array_name = "curr_state",
                    subfunction_prefix = f"fn_{i}",
                    overrides = overrides
                ) 

                # add in all subfunction lines:
                vhdl_str += ("\n    ".join(chaining_lines[:-1]) + "\n    ")

                # add in the MPR feedback
                vhdl_str += (self.component_feedback[i].merge_redundant().generate_VHDL(
                    output_name = f"next_state({i})",
                    array_name = "curr_state",
                    subfunction_prefix = f"fn_{i}",
                    overrides = overrides
                )[0][:-1]) + " XOR "

                # add in the chaining logic
                vhdl_str += chaining_lines[-1][len(f"next_state({i}) <= "):] + "\n    "

            elif self.has_chaining[i] and not include_mpr:
                vhdl_str += ("\n    ".join(self.chaining_feedback[i].merge_redundant().generate_VHDL(
                    output_name = f"next_state({i})",
                    array_name = "curr_state",
                    subfunction_prefix = f"fn_{i}",
                    overrides = overrides
                )) + "\n    ")
                
            elif not self.has_chaining[i] and include_mpr:
                vhdl_str += (self.component_feedback[i].merge_redundant().generate_VHDL(
                    output_name = f"next_state({i})",
                    array_name = "curr_state",
                    subfunction_prefix = f"fn_{i}",
                    overrides = overrides
                )[0] + "\n    ")

            else:
                vhdl_str += (f"next_state({i}) <= '0';\n    ")


            # add in override strings:
            for j, node in enumerate(self.fn_list[i].subfunctions()):
                if node not in overrides:
                    overrides[node] = f'fn_{i}_{j+1}'

        subfunction_vars = list(overrides.values())
        if subfunction_vars:
            subfn_var_string = f"signal {', '.join(overrides.values())}: std_logic;\n"
        else:
            subfn_var_string = ""

        vhdl_str = f"""
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

    signal curr_state, next_state:std_logic_vector({self.size - 1} downto 0);
    {subfn_var_string}    
begin

    statereg: process(i_clk, i_rst)
    begin
        if (i_rst = '1') then
            curr_state <= i_seed_data;
        elsif (i_clk = '1' and i_clk'event) then
            curr_state <= next_state;
        end if;
    end process;\n""" + vhdl_str

        vhdl_str += """
    output <= currstate;

end run;

"""
        with open(filename, "w") as f:
            f.write(vhdl_str)