from ProductRegisters.ANF import ANF

#aaaa broken
from copy import deepcopy
import json

class FeedbackFunction:
    def __init__(self, anf_list):
        #convert update to a list of ANF<int> objects:
        anf_list = [ANF(bitFn for bitFn in anf_list)]
        self.anf = anf_list
        self.size = len(anf_list)

    def __getitem__(self, idx): return self.anf[idx]

    def __setitem__(self, idx, val): self.anf[idx] = val

    def __len__(self): return self.size

    def __eq__(self, other): return self.anf == other.anf

    def __str__(self):
        outstr = ""
        for i in range(len(self.anf)-1,-1,-1):
            outstr += str(i) + "="
            outstr += str(self.anf[i]) + ";\n"
        return outstr[:-1]

    #FUNCTION MANIPULATION:
    #new function that generates the same states in reverse bit order
    #NEEDS UPDATING TO ANF OBJECTS:
    def flip(self):
        newFn = []
        for bitFn in self.anf[::-1]:
            newBitFn = []
            for term in bitFn:
                newTerm = []
                for var in term:
                    newTerm.append(self.size-var-1)
                newBitFn.append(newTerm)
            newFn.append(ANF(newBitFn))
        self.anf = newFn

    """
    #TOO SLOW TO USE
    def compose(self,other):
        from itertools import product, chain
        #other function -> then this function.

        #create appropriate temp/new variables.
        newFunc = []
        for bitFn in self.fn:
            newEq = []
            for term in bitFn:
                #print(f"term: {term}")

                #Handle true terms separately.
                if term == True:
                    newTerms = [True]

                
                else:
                    #gather equations, and replace True -> [] for itertools
                    syms = []
                    for var in term:
                        #print(f"var: {var}; {other.fn[var]}")
                        if other.fn[var]:
                            syms.append([x if not (x is True) else [] for x in other.fn[var]])

                    #print(f"syms: {syms}")
                    
                    #generate new terms from those equations
                    newTerms = (list(chain(*x)) for x in product(*syms))
                    newTerms = (x if x else True for x in newTerms)

                    #print(f"newTerms: {newTerms}")

                #add new terms into the new Eq:
                for nTerm in newTerms:
                    if nTerm in newEq:
                        newEq.remove(nTerm)
                    else:
                        newEq.append(nTerm)
            #print(f"NEWEQ: {newEq}\n")
            newFunc.append(newEq)
        return newFunc
    """     

    #return the number of gates before any optimization (rough estimate of size)
    def gateSummary(self):
        xors = 0
        ands = 0
        for b in self.anf:
            for term in b:
                xors += 1
                ands += len(term)-1
        #print(f"\nMaximum # of XORS: {xors}")
        #print(f"Maximum # of ANDS: {ands}")
        return (xors, ands)

    #return true if the function is linear
    def isLinear(self, allowAfine = False):
        valid = True
        for ls in self.anf:
            for l in ls:
                if type(l) == bool:
                    valid &= allowAfine
                elif len(l) != 1:
                    valid = False
        return valid

    #HELLA BROKEN
    #to/from JSON methods (inherited by all bitfunction types)
    def toJSON(self, filename):
        dct = {}
        dct["Type"] = type(self).__name__
        dct["Data"] = self.__dict__
        with open(filename, "w") as f:
            f.write(json.dumps(dct))

    @classmethod
    def fromJSON(self, filename):
        with open(filename, "r") as f:
            dct = json.load(f)

            #if wrong fn type
            if dct["Type"] != self.__name__:
                expected = self.__name__
                actual = dct["Type"]
                raise ValueError(
                    f"Expected function type {expected}, but got {actual}"
                    )
            
            #otherwise remake the FN
            fn = object.__new__(self)
            fn.__dict__ = dct["Data"]
        return fn

    #writes a VHDL file
    #Credit: Anna Hemingway
    def generateVHDL(self, filename):
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
                writestr = ""
                writestr += f"    nextstate({str(i)}) <= "
                for term in self.anf[i]:
                    if type(term) == bool:
                        writestr += f"currstate({str(i)}) XOR "
                    elif len(term) == 1:
                        writestr += f"currstate({str(term[0])}) XOR "
                    else :
                        writestr += "(" + " AND ".join(f"currstate({str(t)})" for t in term) +") XOR "
                writestr = writestr[:-5] + ";\n"
                f.write(writestr)
            f.write("""

                    output <= currstate;

                end run;

                """)
    
