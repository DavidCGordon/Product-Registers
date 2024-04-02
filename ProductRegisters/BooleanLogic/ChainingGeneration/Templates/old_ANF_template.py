import ProductRegisters.BooleanLogic.ChainingGeneration.TemplateBuilding as tb
from ProductRegisters.BooleanLogic.Gates import *
from ProductRegisters.BooleanLogic.Inputs import *

def old_ANF_function(reference_block, max_and = 4, max_xor = 4):
  n = len(reference_block)

  def first_layer_AND(num_inputs):
    if num_inputs < n:
      k = num_inputs
    else:
      k = n-1

    return (
      tb.DISTINCT(
        parameters = {"group_id": 0, "disable_on_failure": True},
        source = (
          tb.NONCONSTANT(
            parameters = {},
            source = (
              tb.GATE(
                parameters={"gate_class": AND},
                sources = [
                  tb.SAMPLE(
                    parameters = {"num_sampled":k},
                    source_distribution = [
                      (1/n, tb.VALUE(VAR(bit))) for bit in reference_block
                    ]
                  )
                ]
              )
            )
          ) 
        )
      )
    )
  
  def second_layer_XOR(num_inputs):
    return (
      tb.UNIQUE(
        parameters = {"group_id": 0, "disable_on_failure": False},
        source = (
          tb.NONCONSTANT(
            parameters = {},
            source = (
              tb.GATE(
                parameters={"gate_class": XOR},
                sources = [
                  tb.REPEAT(
                    parameters = {"iterations": num_inputs},
                    source = (
                      tb.SAMPLE(
                        parameters = {"num_sampled": 1},
                        source_distribution = [
                          (1/(max_and-1), first_layer_AND(num_inputs = i)) for i in range(2,max_and+1)
                        ]
                      )
                    )
                  )
                ]
              )
            )
          ) 
        )
      )
    )

  n = len(reference_block)
  return tb.SAMPLE(
    parameters = {"num_sampled":1},
    source_distribution = [
      (1/(max_xor-1), second_layer_XOR(num_inputs = i)) for i in range(2,max_xor+1)
    ]
  )


def old_ANF_template(max_and = 4, max_xor = 4):
  def template_fn(cmpr):
    fns = {}
    for i in range(1,cmpr.num_components):
      #template = arman_template(reference_block=C.blocks[i-1])
      template = old_ANF_function(
        reference_block=cmpr.blocks[i-1],
        max_and=max_and, 
        max_xor=max_xor
      )

      for bit in cmpr.blocks[i][::-1]:
        fn = template.sample()[0]
        fns[bit] = fn

    # cleanup
    tb.UNIQUE._function_cache = {}
    tb.DISTINCT._function_cache = {}
    return fns
  return template_fn