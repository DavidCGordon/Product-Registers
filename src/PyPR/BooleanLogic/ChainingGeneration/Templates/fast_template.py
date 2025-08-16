import PyPR.BooleanLogic.ChainingGeneration.TemplateBuilding as tb
from PyPR.BooleanLogic.Gates import *
from PyPR.BooleanLogic.FunctionInputs import *

def fast_function(reference_block, max_and):
  n = len(reference_block)
  if max_and < n:
    k = max_and
  else:
    k = n-1

  return (
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


def fast_template(max_and = 4):
  def template_fn(cmpr):
    fns = {}
    for i in range(1,cmpr.num_components):
      # don't add chaining to TFunction-like segments
      if len(cmpr.blocks[i-1]) == 1 and len(cmpr.blocks[i]) == 1:
        continue
      
      template = fast_function(
        reference_block=set().union(*cmpr.blocks[:i]),
        max_and = max_and
      )

      for bit in cmpr.blocks[i][::-1]:
        fn = template.sample()[0]
        fns[bit] = fn

    return fns
  return template_fn