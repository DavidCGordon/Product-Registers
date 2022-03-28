from itertools import product

# A container class which can hold an ANF of any hashable type
# Note that this class is unordered, because it uses sets. 
class ANF:

    #converts iterables/bools to the internal format of ANF
    def _convert_term(self, term):
        if (type(term) == bool and term == True) or (not term and type(term) != bool):
            return frozenset()
        else:
            return frozenset(term)


    def __init__(self, formula = None):
        #set of frozensets of <basis, count> pairs
        self._xorSet = set()

        #if they want an empty ANF object.
        if formula == None: return

        #convert any nested iterable to proper set/frozenset format:
        for term in formula:

            #accept True, None, (), [], {} for the "1" element, do not accept False.
            self._xorSet = self._xorSet.symmetric_difference([self._convert_term(term)])
        
        #self.varSet = set.union(*(set(tup) for tup in self._xorSet))


    # ANF Operations:

    # ADD and XOR
    def __add__(self,other):return ANF(self._xorSet ^ other._xorSet)
    def __xor__(self,other): return self.__add__(other)

    # MUL and AND 
    # (note that these grow large and shouldn't be overused)
    def __and__(self,other): return self.__mul__(other)
    def __mul__(self,other):
        newANF = ANF()
        for a,b in product(self._xorSet,other._xorSet):
            newANF += ANF([a | b])
        return newANF

    # pretty print for __str__, generic object for repr:
    def __str__(self):
        stringBeginning = ""
        termStrings = []
        for term in self._xorSet:
            if term == frozenset():
                stringBeginning = "True,"
            elif len(term) == 1:
                termStrings.append(repr(tuple(term)[0]) + ",")
            else:

                #if the internal objects have an order, sort terms.
                try: 
                    termStrings.append("(" + ",".join(repr(t) for t in sorted(term)) + "),")
                except TypeError:
                    termStrings.append("(" + ",".join(repr(t) for t in term) + "),")
        
        #sort string terms, and print vaguely by size 
        return stringBeginning + "".join(sorted(termStrings, key = lambda x: len(x)))[:-1]

    # Mimic the set interface.
    def __len__(self): return len(self._xorSet)
    def __eq__(self,other): return self._xorSet == other._xorSet
    def __hash__(self): return self._xorSet.__hash__()
    def __iter__(self): return iter(self._xorSet)
    def add(self,term): self._xorSet.add(self._convert_term(term))
    def remove(self,term): self._xorSet.remove(self._convert_term(term))
    def __contains__(self,term): self._xorSet.__contains__(self._convert_term(term))

