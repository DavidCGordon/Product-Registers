# PyPR — Product Registers

PyPR is a Python research library for simulating and cryptanalyzing feedback
registers: LFSRs (Fibonacci and Galois), NLFSRs, FCSRs, Mersenne product
registers (MPRs), and composite MPRs (CMPRs) built by chaining MPRs together
with nonlinear feedback.

It provides register simulation (compiled with numba or in plain Python),
Boolean function algebra with code generation, period and linear-complexity
analysis, register synthesis from sequences, and algebraic and cube attacks.

## Setup

PyPR needs Python 3.12 or later. From a clone of this repository:

```bash
pip install -e .            # the library
pip install -e ".[dev]"     # plus pytest
pytest                      # run the test suite
```

## Documentation

- [Getting Started](docs/Getting%20Started.md): worked examples.
- [`docs/`](docs/): the theory, architecture and conventions behind the code, indexed by [`docs/Theory_Map.md`](docs/Theory_Map.md).
- [`CLAUDE.md`](CLAUDE.md): contributor conventions.

To browse the docs and an API reference generated from the docstrings as a site:

```bash
pip install -e ".[docs]"
sphinx-build -b html docs docs/_build/html     # then open docs/_build/html/index.html
```

## License

MIT. See [`LICENSE.md`](LICENSE.md).
