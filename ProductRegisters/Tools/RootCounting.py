from itertools import product,chain,combinations
from ProductRegisters.ANF import ANF

#small combinatorial helper fns:
def powerset(iterable):
	    s = list(iterable)
	    return chain.from_iterable(combinations(s, r) for r in range(1,len(s)+1))

def choose(n,k):
    prod = 1
    for i in range(k):
        prod *= (n-i)/(k-i)
    return round(prod)

def binsum(n,k):
    tot = 0
    for i in range(1,k+1):
        tot += choose(n,i)
    return tot

#A small pair class, to designate pairs describing root sets 
class RootPair:
    def __init__(self,basis,count):
        self.basis = basis
        self.count = count

    def __repr__(self):
        return f"<{self.basis},{self.count}>"

#A class which implements general sequence manipulation
class RootExpression:

    #An easy constructor of a given size. 
    @classmethod
    def fromPrimitive(self,n):
        return RootExpression([[RootPair(n,1)]])
    
    
    def __init__(self,rootANF):
        #{tuple: ANF} table.
        #the Label tuple is sorted descending (e.g. [7,5,3]).
        #the Corresponding ANF consists of all terms with the same basis.
        self.rootANF = ANF(rootANF)

    def clone(self):
        newterms = []
        for term in self.rootANF:
            newterms.append([RootPair(rp.basis, rp.count) for rp in term])
        return RootExpression(newterms)

    #construct basis table:
    def basis_table(self):
        table = {}
        for term in self.rootANF:
            label = tuple(sorted([rp.basis for rp in term]))
            if label in table:
                table[label].add(term)
            else:
                table[label] = ANF([term])
        return table 

    #combinationANF :: ANF of ints
    #varMap :: int -> RootExpression map
    @classmethod 
    def fromCombination(self, varMap, combinationANF):
        # initial ANF product
        # this handles "forced" degeneracies. 
        initialANF = ANF()
        for term in combinationANF:
            newTerms = ANF([True])
            for var in term:
                newTerms *= varMap[var].rootANF
            initialANF += newTerms

        # merge same basis pairs and convert full root sets
        # this handles combination logic and overflow.
        convertedANF = ANF()
        for term in initialANF:
            # term is a frozenSet of rootPairs:

            # newTerms holds the ANF of the stuff being added.
            # because it is a product, it is initialized to "1"
            newTerms = ANF([True])

            # group rootPairs by basis (within a single term)
            # table looks like {(basis):[list of RootPairs]}
            table = {}
            for rootpair in term:
                if rootpair.basis in table:
                    table[rootpair.basis] += [rootpair]
                else:
                    table[rootpair.basis] = [rootpair]

            # create newTerms
            for b,ls in table.items():
                # if its a single term, keep the same RootPair for xor-cancellation
                if len(ls) == 1:
                    newTerms *= ANF([[ls[0]]])
                # else, merge counts, and create new terms
                else:
                    c = sum(rp.count for rp in ls)
                    if c >= b:
                        newPair = RootPair(b,b-1)
                        newTerms *= ANF([[newPair],True])
                    else:
                        newPair = RootPair(b,c)
                        newTerms *= ANF([[newPair]])
            convertedANF += newTerms

        # simplify ANF assuming nondegeneracy 
        # NOTE: this introduces error, which grows the more this function is called
        # but, it is also necessary for this to not explode.

        finalANF = ANF()
        table = {}

        #construct basis table using tuples
        for term in convertedANF:

            #create sorted label/count tuples for the term:
            label = []
            counts = []
            ls = sorted([(rp.basis,rp.count) for rp in term], key = lambda x: x[0])
            for basis,count in ls:
                label.append(basis)
                counts.append(count)
            label = tuple(label)
            counts = tuple(counts)

            # if the label is in the table:
            if label in table:

                useful = True
                for other_counts, other_term in table[label]:

                    # if this term's root set is a subset of other term, don't add.
                    if not any((a > b) for a, b in zip(counts,other_counts)):
                        useful = False
                        break

                    # if other root set is a subset of this term's, remove the old one.
                    elif all((a >= b) for a, b in zip(counts,other_counts)):
                        table[label].remove((other_counts,other_term)) 

                #only append to table if it's useful
                if useful:
                    table[label] += [(counts,term)]
            else:
                table[label] = [(counts,term)]

        # convert table back into a single ANF
        for group in table.values():
            finalANF += ANF([term for counts,term in group])
        return RootExpression(finalANF)

    def __add__(self,other):
        return RootExpression(self.rootANF + other.rootANF)
              
    def __str__(self):
        outStr = ""
        for label, anf in self.basis_table().items():
            outStr += "--<"
            empty = True
            for basis in label:
                empty = False
                outStr += f"{basis}, "
            if not empty:
                outStr = outStr[:-2] + ">--\n"
            else:
                outStr += " Identity >--"
        
            for term in anf:
                outStr += "("
                for rp in term:
                    outStr += f"{rp.__str__()}, "
                outStr = outStr[:-2] + ")\n"
        return outStr

    def eval(self):
        LinearComplexity = 0

        #remove duplicates before calculating
        table = {}
        for term in self.rootANF:

            #create sorted label/count tuples for the term:
            label = []
            counts = []
            ls = sorted([(rp.basis,rp.count) for rp in term], key = lambda x: x[0])
            for basis,count in ls:
                label.append(basis)
                counts.append(count)
            label = tuple(label)
            counts = tuple(counts)

            # if the label is in the table()
            if label in table:

                useful = True
                for other_counts in table[label]:
                    #if its less than or equal something already added, break and dont add.
                    if (counts == other_counts) or not any((a > b) for a, b in zip(counts,other_counts)):
                        useful = False
                        break

                    # if its strictly greater than something already added, remove the lesser.
                    # this may not ever proc?
                    elif all((a >= b) for a, b in zip(counts,other_counts)):
                        table[label].remove(other_counts)

                #only append to table if it's useful/unique
                if useful:
                    table[label] += [counts]
            else:
                table[label] = [counts]

        #use relevant terms to calculate LC:             
        for label,group in table.items():
            for subset in powerset(group):
                #calculate intersection:
                intersection = [min(cs) for cs in zip(*subset)]

                #calculate intersection binsum product
                termProduct = 1
                for basis,count in zip(label,intersection):
                    termProduct *= binsum(basis,count)
                sign = (-1)**(len(subset)+1)
                LinearComplexity += (sign*termProduct)
        return LinearComplexity



    """
    def __mul__(self,other):
        # initial ANF product
        initialANF = self.rootANF * other.rootANF

        # merge same basis pairs and
        # convert full root sets
        convertedANF = ANF()
        for term in initialANF:
            # term is a frozenSet of rootPairs:

            # newTerms holds the ANF of the stuff being added.
            # because it is a product, it is initialized to "1"
            newTerms = ANF([True])

            # group rootPairs by basis (within a single term)
            # table looks like {(basis):[list of RootPairs]}
            table = {}
            for rootpair in term:
                if rootpair.basis in table:
                    table[rootpair.basis] += [rootpair]
                else:
                    table[rootpair.basis] = [rootpair]

            # create newTerms
            for b,ls in table.items():
                # if its a single term, keep the same RootPair for xor-cancellation
                if len(ls) == 1:
                    newTerms *= ANF([[ls[0]]])
                # else, merge counts, and create new terms
                else:
                    c = sum(rp.count for rp in ls)
                    if c >= b:
                        newPair = RootPair(b,b-1)
                        newTerms *= ANF([[newPair],True])
                    else:
                        newPair = RootPair(b,c)
                        newTerms *= ANF([[newPair]])
            convertedANF += newTerms
        #remove "useless" terms:
        #construct basis table using tuples
        table = {}
        for term in convertedANF:

            #create sorted label/count tuples for the term:
            label = []
            counts = []
            ls = sorted([(rp.basis,rp.count) for rp in term], key = lambda x: x[0])
            for basis,count in ls:
                label.append(basis)
                counts.append(count)
            label = tuple(label)
            counts = tuple(counts)

            # if the label is in the table:
            if label in table:

                useful = True
                for other_counts, other_term in table[label]:

                    # something in the table is useful according to what we've seen so far
                    #  => if you find something that matches it: stop thinking - just add :)
                    if counts == other_counts:
                        break

                    # if other root set is a subset of this term's, remove the old one.
                    elif all((a >= b) for a, b in zip(counts,other_counts)):
                        table[label].remove((other_counts, other_term))

                    # if this term's root set is a subset of other term.
                    elif not any((a <= b) for a, b in zip(counts,other_counts)) and counts != other_counts:
                        useful = False
                        break

                #only append to table if it's useful
                if useful:
                    table[label] += [(counts,term)]
            else:
                table[label] = [(counts,term)]

        # convert table back into a single ANF
        finalANF = ANF()
        for group in table.values():
            finalANF += ANF([term for counts,term in group])
        return RootExpression(finalANF)
        """
    """
        #create table of identical label/count pairs
        #construct basis table: (basis,count) -> [terms]
        duplicates_table = {}
        for term in convertedANF:
            #create sorted label/count tuples for the term:
            label = []
            counts = []
            ls = sorted([(rp.basis,rp.count) for rp in term], key = lambda x: x[0])
            for basis,count in ls:
                label.append(basis)
                counts.append(count)
            label = tuple(label)
            counts = tuple(counts)
            
            if (label,counts) in duplicates_table:
                duplicates_table[label,counts] += [term]
            else:
                duplicates_table[label,counts] = [term]
        for (label,counts),terms in duplicates_table.items():
            print(label,counts,terms)
            if len(terms) > 1:
                for term in terms:
                    convertedANF.remove(term)
                convertedANF.add(frozenset(RootPair(basis, count) for basis, count in zip(label,counts)))
        
        #return RootExpression(convertedANF)

        #clean ANF of (same basis-lower count) terms
        #create new ANFs grouped by term-basis
        table = {}
        for term in convertedANF:
            label = []
            counts = []
            ls = sorted([(rp.basis,rp.count) for rp in term], key = lambda x: x[0])
            for b,c in ls:
                label.append(b)
                counts.append(c)
            label = tuple(label)
            counts = tuple(counts)
            if label in table:
                table[label] += [(counts,term)]
            else:
                table[label] = [(counts,term)]
        
        #remove duplicates
        for label,group in table.items():
            something_removed = True
            while len(group) > 1 and something_removed:
                for counts1, term1 in group:
                    for counts2, term2 in group:
                        something_removed = False
                        if all((a >= b) for a, b in zip(counts1,counts2)) and term1 != term2:
                            group.remove((counts2,term2))
                            something_removed = True

        #convert lists into ANFs:
        finalANF = ANF()
        for label,group in table.items():
            finalANF += ANF(term for count,term in group)
        """
        #return RootExpression(finalANF
