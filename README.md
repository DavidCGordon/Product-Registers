# Product Registers
This is a set of files designed to build a testing environment for aspects of research into product registers, nlfsrs, and other large period iterated boolean functions.

## Files
BitFunction.py: This file contains the BitFunction class, which is designed to interact with the NLFSR class. 
The BitFunction may be printed to show an ANF representation of each bit's update function with grouped terms being products (AND) and comma separated terms being a sum (XOR).

ProductRegister.py: This file contains the ProductRegister class, which is constructed with a seed value and a BitFunction as its parameters. It is designed to iterate and evalute its given BitFunction over its BitVector contents

StandardMIPR.py: This file contains the code needed to combine an input with a product register iterator to create a multiple input product register (MIPR). this MIPR contains the runHash function, which can be used to test how a product register would evaluate a specific input.  

Various Bitfunction subclasses:
represent different constructions or families of bitfunctions, each of which has its own parameters, and can be made into it's own product registers. many Bitfunctions are constructed using [hexCodes](http://users.ece.cmu.edu/~koopman/lfsr/index.html) associated with primitive polynomials in GF(2), and are expected to be used with the koopman hexcode format.

tools:
the FunctionWriter.py file contains methods for converting any iterated boolean output sequence back into it's ANF format, using the algorithm outlined in [this paper](http://www.selmer.uib.no/odbf/help/ttanf.pdf) and some scaffold code. 

the Tests.py file contains tests for period, as well as an implementation of the berlekamp-massey algorithm and a few other minor useful tests for analyzing product register output. 

## Usage
There are a variety of specialized BitFunction subclasses, which represent different ways to construct useful bitfunctions and have different parameters based on their construction. These can be used to construct Product Registers and MIPRs for testing:

```python
example = ProductRegister(seed, BitFunction(params...))
```

An example with real values, if you want to run for yourself:

```python
example = ProductRegister(seed, CrossJoin(16, "810A")) 
```

## Credit for CrossJoin bitfunctions:
The CrossJoin bitfunction implements the concepts outlined in Elena Dubrova's [A Scalable Method for Constructing Galois NLFSRs with Period 2 <sup>n</sup> −1 using Cross-Join Pairs](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6290394) and 
[A Transformation From the Fibonacci to the Galois NLFSRs](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5290281) to create a generator for randomized NLFSRs, with variables for max term length and tap density. The sampling method used by this implementation for randomizing the CrossJoin terms is my work.

Note on Max-term-length and Tap-density:
while term length is arbitrary, small values like 4 and 5 consistently yield better results in terms of both NLSFR complexity and size, as larger values lead to fewer terms with more gates. Larger term lengths also lead to longer calculation times.

Due to the way The NLFSRs are created, Larger density values  (i.e. 60% - 85%, ideally) tend to yield better results, and calculate faster. Density is represented as a decimal float.