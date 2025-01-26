import ProductRegisters.BooleanLogic.ChainingGeneration.TemplateBuilding as tb
from ProductRegisters.BooleanLogic.ChainingGeneration.Templates.old_ANF_template import old_ANF_function
from ProductRegisters.BooleanLogic.Gates import *
from ProductRegisters.BooleanLogic.Inputs import *

import random

def prob_ANF_template(max_and = 4, max_xor = 4, p = 1.00):
  def template_fn(cmpr):
    fns = {}
    for i in range(1,cmpr.num_components):
      # better integration for T-function like segments:
      if len(cmpr.blocks[i-1]) == 1 and len(cmpr.blocks[i]) == 1:
        continue
      elif len(cmpr.blocks[i-1]) == 1:
        fns[cmpr.blocks[i][-1]] = VAR(cmpr.blocks[i][-1] + 1)
        continue


      template = old_ANF_function(
        reference_block=cmpr.blocks[i-1],
        max_and=max_and, 
        max_xor=max_xor
      )

      for bit in cmpr.blocks[i][::-1]:
        if (random.random() <= p):
          fn = template.sample()[0]
          fns[bit] = fn

    # cleanup
    tb.UNIQUE._function_cache = {}
    tb.DISTINCT._function_cache = {}
    return fns
  return template_fn