import ProductRegisters.BooleanLogic.ChainingGeneration.TemplateBuilding as tb
from ProductRegisters.BooleanLogic.Gates import *
from ProductRegisters.BooleanLogic.Inputs import *

def arman_function(reference_block, max_and):
  n = len(reference_block)
  if max_and < n:
    k = max_and
  else:
    k = n-1

  return (
   tb.GATE(
    parameters={"gate_class": XOR},
    sources=[
      tb.SAMPLE(
        parameters={'num_sampled':1},
        source_distribution = [
          (.5, tb.VALUE(CONST(1))),
          (.5, tb.SAMPLE(
              parameters = {'num_sampled':1},
              source_distribution = [
                (1/n, tb.VALUE(VAR(bit))) for bit in reference_block
              ]
            )
          )
        ]
      ),
      tb.SAMPLE(
        parameters = {'num_sampled': 1},
        source_distribution = [
          (1/n, tb.VALUE(VAR(bit))) for bit in reference_block
        ]
      ),
      tb.UNIQUE(
        parameters = {"group_id": 0, "attempt_limit":50, "disable_on_failure": True},
        source = (
          tb.GATE(
            parameters = {"gate_class": AND},
            sources = [
              tb.SAMPLE(
                parameters = {'num_sampled': k},
                source_distribution = [
                  (1/n, tb.VALUE(VAR(bit))) for bit in reference_block
                ]
              )
            ]
          )
        )
      )
    ]
  ))


def arman_template(max_and = 4):
  def template_fn(cmpr):
    fns = {}
    for i in range(1,cmpr.num_components):
      template = arman_function(
        reference_block=cmpr.blocks[i-1],
        max_and = max_and
      )

      for bit in cmpr.blocks[i][::-1]:
        fn = template.sample()[0]
        fns[bit] = fn

    # cleanup
    tb.UNIQUE._function_cache = {}
    tb.DISTINCT._function_cache = {}
    return fns
  return template_fn