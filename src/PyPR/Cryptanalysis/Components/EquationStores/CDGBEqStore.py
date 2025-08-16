from PyPR.Cryptanalysis.Components.EquationStores.GrobnerEqStore import GroebnerEqStore
from PyPR.BooleanLogic import XOR, VAR, CONST

class CDGBEqStore:
    def __init__(self):
        self.s0 = GroebnerEqStore()
        self.s1 = GroebnerEqStore()
        self.assumptions = []
        self.var_counts = {}

    def branch_from(self, idx):
        # determine copy direction
        if idx == 0:
            source = self.s0
            dest = self.s1
        elif idx == 1:
            source = self.s1
            dest = self.s0
        else:
            raise ValueError('invalid idx')

        # copy data
        dest.num_vars = source.num_vars
        dest.num_eqs = source.num_eqs
        dest.seen = source.seen.copy()
        dest.queue = source.queue.copy()
        dest.equations = source.equations.copy() 
        dest.lead_terms = source.lead_terms.copy()
        dest.unknown_vars = source.unknown_vars.copy()
        dest.solved_vars = source.solved_vars.copy()

        # select new variables
        unknown_vars = (self.var_counts.keys() - source.solved_vars.keys())
        if unknown_vars:
            selected_var = max(unknown_vars, key= lambda v: self.var_counts[v])
            self.assumptions = [
                XOR(VAR(selected_var), CONST(0)),
                XOR(VAR(selected_var), CONST(1))
            ]
            self.s0.enqueue_equation(self.assumptions[0])
            self.s1.enqueue_equation(self.assumptions[1]) 


    def enqueue_equation(self, equation, extra_const = 0, identifier=None, translate_ANF = True):
        self.var_counts.update({v:1 for v in equation.idxs_used()})
        self.s0.enqueue_equation(equation,extra_const,identifier,translate_ANF)
        self.s1.enqueue_equation(equation,extra_const,identifier,translate_ANF)
        
    def consume_queue(self, num=None, verbose = False):
        if num == None:
            num = -1
        
        while num:
            # set assumptions if this is the first equation:
            if not self.assumptions and len(self.s0.equations) >= 100:
                print("INIT")
                self.branch_from(0)

            # update counts:
            for v in self.s0.queue[0].lead_term:
                self.var_counts[v] += 1
            for v in self.s1.queue[0].lead_term:
                self.var_counts[v] += 1

            # attempt to consume:
            if not (self.s0.consume_queue(1,verbose)):
                print("Branch from 0")
                self.branch_from(1)
            if not (self.s1.consume_queue(1,verbose)):
                print("Branch from 1")
                self.branch_from(0)
            num -= 1
        return True
    
    @property
    def solved_vars(self): 
        return dict(self.s0.solved_vars.items() & self.s1.solved_vars.items())
    
    @property
    def potential_solution(self):
        return (
            (self.s0.solved_vars.keys() == self.var_counts.keys()) or 
            (self.s1.solved_vars.keys() == self.var_counts.keys())
        )