import ProductRegisters.BooleanLogic.ChainingGeneration.TemplateBuilding as tb
from ProductRegisters.BooleanLogic.Gates import *
from ProductRegisters.BooleanLogic.Inputs import *

MAJ3 = (
  XOR(
    AND(VAR(0),VAR(1)),
    AND(VAR(0),VAR(2)),
    AND(VAR(1),VAR(2))
  )
)


def maj_function(reference_block,algebraic_degree,correlation_immunity,require_unique = False):
  n = len(reference_block)
  def xor_leaf():
    return (
      tb.GATE(
        parameters={"gate_class":XOR},
        sources = [
          tb.REPEAT(
            parameters={"iterations": correlation_immunity+1},
            source=(
              tb.DISTINCT(
                parameters={"group_id":0},
                source = (
                  tb.SAMPLE(
                    parameters = {"num_sampled":1},
                    source_distribution = [
                      (1/n, tb.VALUE(VAR(bit))) for bit in reference_block
                    ]
                  )
                )
              )
            )
          )
        ]
      )
    )
  
  def maj_node(degree_remaining):
    # BASE CASE:
    if degree_remaining == 1:
      return xor_leaf()

    # RECURSIVE CASE:
    # branch 1 = largest power of 2 stricly less than the remaining
    branch_1 = 1
    while branch_1 < degree_remaining:
      branch_1 *= 2
    branch_1 //= 2

    # branch 2 = whatever is left
    branch_2 = degree_remaining - branch_1

    return (
      tb.FUNCTION(
        parameters={"fn": MAJ3},
        sources= [
          maj_node(branch_1),
          maj_node(branch_2),
          xor_leaf()
        ]
      )
    )

  if require_unique:
    return (
      tb.UNIQUE(
        parameters = {"group_id": 1},
        source = maj_node(algebraic_degree)
      )
    )
  else:
    return maj_node(algebraic_degree)


def num_necessary_vars(algebraic_degree,correlation_immunity):
  return (correlation_immunity + 1)*(2*algebraic_degree-1)

# default value only works when chaining from registers at least 31
def three_majority_template(correlation_immunity = 3, algebraic_degree = 4, require_unique= False):
  def template_fn(cmpr):
    # input handling:
    nonlocal correlation_immunity
    nonlocal algebraic_degree

    if type(correlation_immunity) == int:
      correlation_immunity = [correlation_immunity]*cmpr.num_components
    if type(algebraic_degree) == int:
      algebraic_degree = [algebraic_degree]*cmpr.num_components
    
    for block_idx in range(cmpr.num_components-1):
      threshold = num_necessary_vars(algebraic_degree[block_idx],correlation_immunity[block_idx])
      if len(cmpr.blocks[block_idx]) < threshold:
        raise ValueError(
          "Unable to complete 3-Maj chaining generation: " +
          f"Requested correlation immunity {correlation_immunity[block_idx]} and "
          f"algebraic degree {algebraic_degree[block_idx]}, which is not possible with "
          f"block size {len(cmpr.blocks[block_idx])} ({threshold} required)."
        )

    # actually making the functions:
    fns = {}
    for block_idx in range(cmpr.num_components-1):
      # better integration for T-function like segments
      if len(cmpr.blocks[block_idx-1]) == 1 and len(cmpr.blocks[block_idx]) == 1:
        continue
      elif len(cmpr.blocks[block_idx-1]) == 1:
        fns[cmpr.blocks[block_idx][-1]] = VAR(cmpr.blocks[block_idx][-1] + 1)


      template = maj_function(
        reference_block=cmpr.blocks[block_idx],
        algebraic_degree=algebraic_degree[block_idx], 
        correlation_immunity=correlation_immunity[block_idx],
        require_unique=require_unique
      )

      for bit in cmpr.blocks[block_idx+1][::-1]:
        fn = template.sample()[0]
        fns[bit] = fn

    # cleanup
    tb.UNIQUE._function_cache = {}
    tb.DISTINCT._function_cache = {}
    return fns
  return template_fn