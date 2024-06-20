from ProductRegisters.Tools.RootCounting.Combinatorics import binsum

def optimistic_evaluation(b,c):
    b,c = int(b), int(c)
    return binsum(b,c)

def pessimistic_evaluation(b,c):
    b,c = int(b), int(c)
    if b == c: 
        return binsum(b,c-1)
    else:
        return binsum(b,c)

def pessimistic_expected_value(b,c):
    b,c = int(b), int(c)
    if b == c: 
        return (binsum(b,c-1) * (2**b-1)**2) / (2**(2*b))
    else:
        return (binsum(b,c) * (2**b-1)**2) / (2**(2*b))

""" 
WORK IN PROGRESS DO NOT USE:
def confidence_interval_lower(b,c,base_function, z = 2.58):
    mean = base_function(b,c) * (2**b-1)**2) / (2**(2*b))
    deviation = sqrt(((2**b-1) / 2**(2*b)) / (base_function(b,c) / b))
    return mean - z*deviation

def confidence_interval_upper(b,c,base_function, z = 2.58):
    mean = base_function(b,c) * (2**b-1)**2) / (2**(2*b))
    deviation = sqrt(((2**b-1) / 2**(2*b)) / (base_function(b,c) / b))
    return mean + z*deviation
"""