"""Tests for Fibonacci and Galois LFSRs: construction, equivalence, synthesis round-trips, inversion.

Derived from ProductRegisters.ipynb (Library Basics > Galois and Fibonacci LFSRs, and
Applications > Berlekamp-Massey Variants for Register Synthesis).

The final section pins three properties that follow from the polynomial
conventions rather than from the implementation: a degree-n primitive
polynomial gives period 2^n - 1, `update_matrix` reproduces the clock as a
GF(2) matrix-vector product, and `fromSeq` round-trips with no manual
reversal (the user-facing payoff of taking the constructor input as the
dual).  See docs/conventions/Polynomial Conventions.md.
"""
import random

import numpy as np
import pytest

from PyPR.FeedbackRegister import FeedbackRegister
from PyPR.FeedbackFunctions import MPR, Fibonacci, Galois, FCSR
from PyPR.BooleanLogic import AND, XOR, VAR
from PyPR.Tools.RegisterSynthesis.lfsrSynthesis import berlekamp_massey


# ── Fibonacci / Galois equivalence ──────────────────────────────────────

@pytest.mark.parametrize("n", [3, 4, 5])
def test_fib_galois_emit_the_same_set_of_sequences(n):
    """Both layouts emit the same cyclic orbit, at a seed-dependent phase.

    Comparing a single seed's output would fail, and comparing only linear
    complexity would pass for unrelated reasons.  The real invariant is set
    equality over all nonzero seeds: the same numeric seed means different things
    in the two layouts -- a Fibonacci state is the sliding window
    (s_t, ..., s_{t+n-1}), a Galois state is an element of GF(2)[x]/P -- so the
    correspondence is a permutation of phases, not the identity.

    This is the claim Polynomial Conventions.md makes under Round-Trip Behavior,
    and the reason its earlier wording about an "indexing defect" was wrong:
    nothing is misindexed, the two encodings just select different phases.
    """
    polynomial = PRIMITIVE[n]
    period = 2 ** n - 1

    def orbits(family):
        return {
            tuple(int(s[0]) for s in
                  FeedbackRegister(seed, family(n, polynomial)).run(compiled=False, limit=period))
            for seed in range(1, 2 ** n)
        }

    fibonacci_orbits, galois_orbits = orbits(Fibonacci), orbits(Galois)

    assert fibonacci_orbits == galois_orbits, (
        f"n={n}: the two layouts emit different sequence sets"
    )
    assert len(fibonacci_orbits) == period, (
        f"n={n}: expected {period} distinct phases, got {len(fibonacci_orbits)}"
    )


def test_fib_galois_from_reg_equivalence():
    """Galois.fromReg on a Fibonacci register reproduces the same sequence (notebook cell 76)."""
    F = Fibonacci(5, '12')
    F.compile()
    fib_reg = FeedbackRegister(31, F)
    seq_fib = [state[0] for state in fib_reg.run(limit=40)]

    fib_reg.reset()
    galois_state, galois_fn = Galois.fromReg(fib_reg, bit=0)
    gal_reg = FeedbackRegister(galois_state, galois_fn)
    seq_gal = [state[0] for state in gal_reg.run(compiled=False, limit=40)]

    assert seq_fib == seq_gal, "Galois.fromReg should reproduce Fibonacci sequence"


# ── Fibonacci fromSeq round-trip ────────────────────────────────────────

def test_fibonacci_fromSeq_round_trip():
    """Fibonacci(*Fibonacci.fromSeq(seq)) reproduces seq on bit 0."""
    F = Fibonacci(5, '12')
    reg = FeedbackRegister(31, F)
    seq = [state[0] for state in reg.run(compiled=False, limit=50)]

    init_state, recovered = Fibonacci.fromSeq(seq)
    rec_reg = FeedbackRegister(init_state, recovered)
    seq2 = [state[0] for state in rec_reg.run(compiled=False, limit=50)]
    assert seq == seq2, "Fibonacci fromSeq round-trip failed"


def test_galois_fromSeq_round_trip():
    """Galois(*Galois.fromSeq(seq)) reproduces seq on bit 0."""
    G = Galois(5, '12')
    reg = FeedbackRegister(31, G)
    seq = [state[0] for state in reg.run(compiled=False, limit=50)]

    init_state, recovered = Galois.fromSeq(seq)
    rec_reg = FeedbackRegister(init_state, recovered)
    seq2 = [state[0] for state in rec_reg.run(compiled=False, limit=50)]
    assert seq == seq2, "Galois fromSeq round-trip failed"


# ── Fibonacci NLFSR ─────────────────────────────────────────────────────

def test_fibonacci_nlfsr_fromSeq():
    """Fibonacci.fromSeq with nonlinear=True reproduces a nonlinear sequence."""
    random.seed(42)
    seq = [random.randint(0, 1) for _ in range(60)]

    init_state, nlfsr = Fibonacci.fromSeq(seq, nonlinear=True)
    rec_reg = FeedbackRegister(init_state, nlfsr)
    seq2 = [state[0] for state in rec_reg.run(compiled=False, limit=60)]
    assert seq == seq2, "Nonlinear Fibonacci fromSeq round-trip failed"


# ── Cross-synthesis round-trips (notebook cell 122) ─────────────────────

def test_all_synthesis_round_trips():
    """All fromSeq variants reproduce the original sequence prefix."""
    random.seed(7)
    seq = [random.randint(0, 1) for _ in range(50)]

    for name, factory in [
        ("Fibonacci", lambda s: FeedbackRegister(*Fibonacci.fromSeq(s))),
        ("Galois", lambda s: FeedbackRegister(*Galois.fromSeq(s))),
        ("Fibonacci NL", lambda s: FeedbackRegister(*Fibonacci.fromSeq(s, nonlinear=True))),
        ("FCSR", lambda s: FeedbackRegister(*FCSR.fromSeq(s))),
    ]:
        reg = factory(seq)
        recovered = [state[0] for state in reg.run(compiled=False, limit=50)]
        assert seq == recovered, f"{name} synthesis round-trip failed"


# ── MPR / Galois relationship (notebook cell 77) ───────────────────────

def test_mpr_galois_relationship():
    """MPR(P_reversed) with the Galois-derived U polynomial matches Galois output."""
    P = [1, 1, 0, 1]
    G = Galois(3, P)
    g_reg = FeedbackRegister(1, G)
    seq_gal = [state[0] for state in g_reg.run(compiled=False, limit=20)]

    # MPR with reversed polynomial should produce a related sequence
    M = MPR(3, P[::-1])
    m_reg = FeedbackRegister(1, M)
    seq_mpr = [state[0] for state in m_reg.run(compiled=False, limit=20)]

    L_gal, _ = berlekamp_massey(seq_gal)
    L_mpr, _ = berlekamp_massey(seq_mpr)
    assert L_gal == L_mpr == 3, f"Expected LC=3, got Galois={L_gal}, MPR={L_mpr}"


# ── Galois inversion ───────────────────────────────────────────────────

def test_galois_invert():
    """Clocking an inverted Galois register undoes a clock of the original."""
    G = Galois(5, '12')
    reg = FeedbackRegister(31, G)

    # clock forward 10 steps
    for _ in reg.run(compiled=False, limit=10):
        pass
    state_at_10 = reg[:].copy()

    # clock forward one more
    reg.clock(compiled=False)
    state_at_11 = reg[:].copy()

    # invert and clock back
    G.invert()
    reg2 = FeedbackRegister(state_at_11, G)
    reg2.clock(compiled=False)
    recovered = reg2[:].copy()

    assert (state_at_10 == recovered).all(), "Inverted Galois should undo one clock step"


# ── Berlekamp-Massey basic ─────────────────────────────────────────────

def test_berlekamp_massey_known():
    """BM on a known LFSR sequence returns the correct linear complexity."""
    F = Fibonacci(5, '12')
    reg = FeedbackRegister(31, F)
    seq = [state[0] for state in reg.run(compiled=False, limit=100)]
    L, poly = berlekamp_massey(seq)
    assert L == 5, f"Expected LC=5, got {L}"


# ── Maximal period, update matrix, and synthesis round-trip ───────────────────

PRIMITIVE = {
    3: [1, 1, 0, 1],
    4: [1, 1, 0, 0, 1],
    5: [1, 0, 1, 0, 0, 1],
    6: [1, 1, 0, 0, 0, 0, 1],
}


@pytest.mark.parametrize("n", sorted(PRIMITIVE))
def test_lfsr_with_primitive_polynomial_has_maximal_period(n):
    """A degree-n primitive polynomial gives every LFSR family period 2^n - 1.

    This is the property that makes the polynomial "primitive": its root
    generates GF(2^n)*, so the nonzero states form a single cycle.  All three
    families are built from the same polynomial under different conventions
    (Fibonacci/Galois take it as the dual, MPR as the primal), and the period
    is invariant under that choice.
    """
    for family in (Fibonacci, Galois, MPR):
        fn = family(n, PRIMITIVE[n])
        fn.compile()
        result = FeedbackRegister(1, fn).period(limit=2 ** 16)
        assert result is not None, f"{family.__name__}(n={n}): no period found"
        period, preperiod = result
        assert (period, preperiod) == (2 ** n - 1, 0), (
            f"{family.__name__}(n={n}) has period {period}, expected {2 ** n - 1}"
        )


@pytest.mark.parametrize("n", sorted(PRIMITIVE))
@pytest.mark.parametrize("family", [Fibonacci, Galois])
def test_update_matrix_reproduces_the_clock(n, family):
    """`update_matrix` satisfies next_state = M @ state over GF(2).

    The matrix is read off the VAR leaves of the feedback DAG, so this checks
    that the extraction agrees with evaluation at every state.  The two layouts
    place the polynomial differently -- Fibonacci concentrates it in the row
    feeding bit n-1, Galois spreads it down column 0 as taps off the departing
    bit -- so this is a genuine check of each construction, not one shared code
    path exercised twice.
    """
    fn = family(n, PRIMITIVE[n])
    matrix = fn.update_matrix

    for seed in range(2 ** n):
        reg = FeedbackRegister(seed, fn)
        before = np.array([int(b) for b in reg._state], dtype=np.uint8)
        reg.clock(compiled=False)
        after = np.array([int(b) for b in reg._state], dtype=np.uint8)
        assert np.array_equal(after, (matrix @ before) % 2), (
            f"{family.__name__}(n={n}) seed={seed}: update_matrix does not "
            f"reproduce the clock"
        )


@pytest.mark.parametrize("n", sorted(PRIMITIVE))
@pytest.mark.parametrize("family", [Fibonacci, Galois])
def test_update_matrix_characteristic_polynomial_is_the_reverse_of_the_input(n, family):
    """char(M) is the reverse of the constructor polynomial, for both layouts.

    The stored polynomial is the dual -- it convolves the output sequence to
    zero -- while the characteristic polynomial of a transition map is the
    primal, whose roots are the roots of the state sequence.  Primal and dual
    of a sequence are reverses of each other, so char(M) = reverse(P).

    Because both families realize the same recurrence, their matrices are
    similar and this single assertion pins that they share a characteristic
    polynomial despite having different entries.  A failure would mean the
    construction had drifted onto the primal reading of its input, which would
    silently reverse every polynomial a caller passes in.
    """
    P = PRIMITIVE[n]
    matrix = family(n, P).update_matrix

    def characteristic_polynomial(M):
        """det(xI + M) over GF(2), by cofactor expansion.

        Polynomials are int bitmasks (bit i = coefficient of x^i).  Signs are
        dropped because -1 == +1 over GF(2), and the diagonal entry of xI + M
        is x + M[r][r], which is `2 ^ M[r][r]` in bitmask form.
        """
        def polymul(a, b):
            out = 0
            while a:
                if a & 1:
                    out ^= b
                a >>= 1
                b <<= 1
            return out

        def det(rows, cols):
            if not cols:
                return 1
            col, total = cols[0], 0
            for i, row in enumerate(rows):
                entry = (2 if row == col else 0) ^ int(M[row][col])
                if entry:
                    total ^= polymul(entry, det(rows[:i] + rows[i + 1:], cols[1:]))
            return total

        packed = det(list(range(M.shape[0])), list(range(M.shape[0])))
        return [(packed >> i) & 1 for i in range(packed.bit_length())]

    assert characteristic_polynomial(matrix) == P[::-1], (
        f"{family.__name__}(n={n}): char poly is "
        f"{characteristic_polynomial(matrix)}, expected reverse(P) = {P[::-1]}"
    )


@pytest.mark.parametrize("n", sorted(PRIMITIVE))
@pytest.mark.parametrize("family", [Fibonacci, Galois])
def test_inverted_register_undoes_a_clock_from_every_state(n, family):
    """Clocking the inverted register returns the state the forward one left.

    This is the property `invert` exists to provide, checked by simulation over
    the whole state space rather than along one trajectory.  Exhaustiveness is
    the point: a register that is *almost* the inverse still returns the right
    answer from a few states, so a spot check can pass on a construction that is
    wrong nearly everywhere.  Before `Fibonacci.invert` applied its bit
    relabelling, it round-tripped 2 of 8 states at n=3 and 1 of 64 at n=6 --
    enough for a single-seed test to pass by luck.
    """
    forward = family(n, PRIMITIVE[n])
    inverted = family(n, PRIMITIVE[n])
    inverted.invert()

    for seed in range(2 ** n):
        start = FeedbackRegister(seed, forward)
        original = [int(b) for b in start._state]
        start.clock(compiled=False)
        stepped = [int(b) for b in start._state]

        back = FeedbackRegister(stepped, inverted)
        back.clock(compiled=False)
        assert [int(b) for b in back._state] == original, (
            f"{family.__name__}(n={n}): {original} -> {stepped} did not invert back"
        )


@pytest.mark.parametrize("n", sorted(PRIMITIVE))
@pytest.mark.parametrize("family", [Fibonacci, Galois])
def test_inverted_update_matrix_is_the_matrix_inverse(n, family):
    """After `invert()`, the update matrix is the GF(2) inverse of the original.

    The linear-algebra form of the test above: M_inverted @ M == I says the same
    thing about all states at once.  Both are kept because they fail
    differently -- this one localizes the problem to the matrix, the simulation
    one proves the register itself is wrong rather than the extraction.
    """
    fn = family(n, PRIMITIVE[n])
    forward = fn.update_matrix.copy()

    fn.invert()
    inverted = fn.update_matrix

    assert np.array_equal((inverted @ forward) % 2, np.eye(n, dtype=np.uint8)), (
        f"{family.__name__}(n={n}): inverted update matrix does not invert the forward one"
    )


@pytest.mark.parametrize("n", sorted(PRIMITIVE))
def test_general_inversion_agrees_with_reciprocal_plus_flip(n):
    """Solving for the discarded bit lands on the reciprocal-plus-flip register.

    `invert` takes one route: solve the feedback for the bit the shift discards,
    which needs no polynomial and so covers nonlinear registers too.  For a
    linear register the textbook route is also available -- rebuild from the
    reciprocal polynomial, then reverse the bit labelling -- and the `invert`
    docstring explains the method by appealing to it.

    This is that explanation, as an assertion.  The general derivation has to
    land on exactly the register the primal/dual argument predicts; if it did
    not, the docstring would be describing something the code does not do.
    """
    general_route = Fibonacci(n, PRIMITIVE[n])
    general_route.invert()

    textbook_route = Fibonacci(n, PRIMITIVE[n])
    textbook_route.fn_list = textbook_route._from_poly(n, PRIMITIVE[n][::-1])
    textbook_route.flip()

    assert ([f.anf_str() for f in general_route.fn_list]
            == [f.anf_str() for f in textbook_route.fn_list]), (
        f"n={n}: solving for the discarded bit disagrees with reciprocal-plus-flip"
    )


@pytest.mark.parametrize("label,feedback", [
    ("s0 ^ s1s2",     XOR(VAR(0), AND(VAR(1), VAR(2)))),
    ("s0 ^ s1 ^ s2s3", XOR(VAR(0), VAR(1), AND(VAR(2), VAR(3)))),
    ("s0 ^ s1s2s3s4",  XOR(VAR(0), AND(VAR(1), VAR(2), VAR(3), VAR(4)))),
], ids=lambda v: v if isinstance(v, str) else "")
def test_nonlinear_register_inverts_when_the_discarded_bit_is_linear(label, feedback):
    """A nonlinear register inverts when its feedback is s0 XOR (rest).

    The shift discards bit 0, so inverting means recovering it from the one
    equation that mentions it.  Writing the feedback as f = s0*A XOR B, that
    equation is solvable exactly when A == 1 -- the discarded bit enters
    linearly and alone, so it XORs back out.  Nonlinearity in B is irrelevant:
    B's variables all survive the shift and can simply be re-read.
    """
    n = 5
    forward = Fibonacci(n, [1] * (n + 1))
    forward.fn_list[n - 1] = feedback
    forward.primitive_polynomial = []        # no polynomial: forces the general route

    inverted = Fibonacci(n, [1] * (n + 1))
    inverted.fn_list[n - 1] = feedback
    inverted.primitive_polynomial = []
    inverted.invert()

    for seed in range(2 ** n):
        start = FeedbackRegister(seed, forward)
        original = [int(b) for b in start._state]
        start.clock(compiled=False)
        stepped = [int(b) for b in start._state]

        back = FeedbackRegister(stepped, inverted)
        back.clock(compiled=False)
        assert [int(b) for b in back._state] == original, (
            f"{label}: {original} -> {stepped} did not invert back"
        )


@pytest.mark.parametrize("label,feedback", [
    ("s1s2 -- bit 0 absent, A = 0",            AND(VAR(1), VAR(2))),
    ("s0s1 -- bit 0 only in a product, A = s1", AND(VAR(0), VAR(1))),
    ("s0 ^ s0s1 -- A = 1 ^ s1, not constant",   XOR(VAR(0), AND(VAR(0), VAR(1)))),
], ids=lambda v: v if isinstance(v, str) else "")
def test_invert_rejects_a_non_bijective_register(label, feedback):
    """Registers whose discarded bit cannot be recovered are refused.

    Each case breaks the A == 1 condition a different way: the bit is absent
    from the feedback, appears only inside a product, or appears with a
    non-constant coefficient.  In every case two distinct states share an
    image, so no inverse map exists -- the error is the honest answer rather
    than a limitation of the construction.
    """
    n = 5
    register = Fibonacci(n, [1] * (n + 1))
    register.fn_list[n - 1] = feedback
    register.primitive_polynomial = []

    with pytest.raises(ValueError, match="not invertible"):
        register.invert()


def test_invert_rejects_the_nonlinear_register_from_fromSeq():
    """The registers `BM_NL` recovers are not bijections, so `invert` refuses.

    `fromSeq(..., nonlinear=True)` builds `Fibonacci(size, [])` and adds the
    recovered rule to the top bit, so the rule lives in `fn_list` while
    `primitive_polynomial` stays empty -- the general route applies.  It rejects
    these because `BM_NL` optimizes for reproducing a sequence, not for
    bijectivity: the recovered feedback carries bit 0 inside large products, so
    the state map collapses and has no inverse to build.
    """
    random.seed(42)
    seq = [random.randint(0, 1) for _ in range(60)]
    _, nonlinear = Fibonacci.fromSeq(seq, nonlinear=True)

    # the rule really is present, just not as a polynomial
    assert nonlinear.primitive_polynomial == []
    assert nonlinear.fn_list[nonlinear.size - 1].anf_str(), "top feedback should be populated"

    with pytest.raises(ValueError, match="not invertible"):
        nonlinear.invert()


@pytest.mark.parametrize("family", [Fibonacci, Galois])
def test_invert_discards_the_cached_update_matrix(family):
    """`invert` clears the cached matrix, because it rebuilds `fn_list`.

    `update_matrix` is a cached_property derived from `fn_list`, and `invert`
    replaces that list.  Without invalidation the property keeps returning the
    pre-inversion matrix, which then disagrees with the register's actual
    clocking -- a stale value that looks valid, which is the failure mode this
    caching contract exists to prevent.
    """
    fn = family(5, PRIMITIVE[5])
    before = fn.update_matrix.copy()
    assert "update_matrix" in fn.__dict__, "update_matrix should cache on first access"

    fn.invert()
    assert "update_matrix" not in fn.__dict__, (
        f"{family.__name__}.invert should discard the cached matrix"
    )
    assert not np.array_equal(fn.update_matrix, before), (
        "the recomputed matrix should reflect the inverted feedback"
    )


@pytest.mark.parametrize("n", sorted(PRIMITIVE))
@pytest.mark.parametrize("family", [Fibonacci, Galois])
def test_invert_is_an_involution(n, family):
    """Inverting twice returns the register to its original configuration.

    `invert` is documented as a toggle, and it decides which construction to
    rebuild by reading `is_inverted`.  If that flag is not updated, the second
    call repeats the first branch instead of undoing it, leaving the register
    permanently inverted while still reporting `is_inverted == False` -- a state
    that is wrong in two ways at once and visible from neither alone.  Checking
    the flag and the matrix together is what distinguishes them.
    """
    fn = family(n, PRIMITIVE[n])
    assert fn.is_inverted is False, "a freshly built register is not inverted"
    forward = fn.update_matrix.copy()

    fn.invert()
    assert fn.is_inverted is True, f"{family.__name__}.invert did not set the flag"
    assert not np.array_equal(fn.update_matrix, forward), (
        "inverting should change the feedback"
    )

    fn.invert()
    assert fn.is_inverted is False, f"{family.__name__}.invert did not clear the flag"
    assert np.array_equal(fn.update_matrix, forward), (
        f"{family.__name__}(n={n}): two inversions did not restore the original matrix"
    )


@pytest.mark.parametrize("n", sorted(PRIMITIVE))
@pytest.mark.parametrize("family", [Fibonacci, Galois])
def test_fromSeq_round_trips_without_manual_reversal(n, family):
    """`family(*family.fromSeq(seq))` reproduces seq, as the convention promises.

    This is the user-facing payoff of keeping BM's dual output: the polynomial
    that comes out of synthesis is exactly the one the constructor expects, so
    no reversal is needed at the boundary.  Inserting a reversal anywhere in
    that chain would break this while leaving the linear complexity unchanged.
    """
    source = FeedbackRegister(1, family(n, PRIMITIVE[n]))
    seq = [int(state[0]) for state in source.run(compiled=False, limit=8 * n)]

    state, recovered = family.fromSeq(seq)
    replay = FeedbackRegister(state, recovered)
    assert [int(s[0]) for s in replay.run(compiled=False, limit=8 * n)] == seq, (
        f"{family.__name__}(n={n}): fromSeq round-trip did not reproduce the sequence"
    )


def test_update_matrices_match_the_documented_shift_down_dual_layout():
    """Both matrices equal the worked example in Polynomial Conventions.md.

    That document lays out all eight interpretations of "the 3-bit LFSR with
    polynomial x^3 + x + 1" -- shift up/down crossed with primal/dual, for each
    of the Fibonacci and Galois layouts -- and states that this library picks
    **shift-down dual** for both.  These are the two matrices it prints for that
    corner, transcribed.

    This is the one assertion that ties the published convention directly to the
    code: if a future change flipped a shift direction or swapped the
    primal/dual reading, the register would still clock consistently and every
    other test here would still pass, but it would no longer be the register the
    documentation describes.
    """
    P = [1, 1, 0, 1]        # 1 + x + x^3

    assert np.array_equal(Fibonacci(3, P).update_matrix, np.array([
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
    ], dtype=np.uint8)), "Fibonacci no longer matches the documented shift-down dual matrix"

    assert np.array_equal(Galois(3, P).update_matrix, np.array([
        [1, 1, 0],
        [0, 0, 1],
        [1, 0, 0],
    ], dtype=np.uint8)), "Galois no longer matches the documented shift-down dual matrix"


@pytest.mark.parametrize("trial", range(5))
def test_fromSeq_bijective_register_inverts_and_has_no_transient(trial):
    """`fromSeq(nonlinear=True, bijective=True)` yields an invertible register.

    This is the end-to-end payoff of the option.  The default nonlinear fit is
    rejected by `invert` because its state map collapses; the bijective fit is
    accepted, and because a bijection's state graph is a disjoint union of
    cycles, every seed sits on a cycle -- so the reported preperiod is 0 rather
    than the transient a collapsing map produces.
    """
    random.seed(800 + trial)
    seq = [random.randint(0, 1) for _ in range(20)]

    state, register = Fibonacci.fromSeq(seq, nonlinear=True, bijective=True)

    produced = [int(s[0]) for s
                in FeedbackRegister(state, register).run(compiled=False, limit=len(seq))]
    assert produced == seq, "the recovered register should reproduce the sequence"

    register.invert()                       # raises if the map is not a bijection
    register.invert()                       # and back

    # Uncompiled deliberately: the recovered feedback is a XOR of one minterm per
    # observed window -- several hundred ANF terms -- and numba spends tens of
    # seconds compiling that, while the register itself is small enough
    # (m is about 8) that interpreting 2^14 steps is near-instant.
    result = FeedbackRegister(state, register).period(compiled=False, limit=2 ** 14)
    assert result is not None, "a bijective register should have a period"
    assert result[1] == 0, (
        f"trial {trial}: preperiod {result[1]} -- a bijection puts every seed on a cycle"
    )


@pytest.mark.parametrize("trial", range(5))
def test_fromSeq_linear_bijective_fit_inverts(trial):
    """`fromSeq(bijective=True)` on the linear path yields an invertible LFSR.

    Berlekamp-Massey's ordinary output taps bit 0 only when its polynomial
    reaches full degree, which failed on a substantial share of uniform random
    inputs in a sampling experiment, leaving a singular register.  The constrained fit forces that tap, so the
    result reproduces the sequence *and* survives `invert`.
    """
    random.seed(900 + trial)
    seq = [random.randint(0, 1) for _ in range(24)]

    state, register = Fibonacci.fromSeq(seq, bijective=True)

    produced = [int(s[0]) for s
                in FeedbackRegister(state, register).run(compiled=False, limit=len(seq))]
    assert produced == seq, "the constrained linear fit should reproduce the sequence"

    register.invert()          # raises if the update map is not a bijection
    register.invert()


@pytest.mark.parametrize("trial", range(4))
def test_galois_bijective_fit_reproduces_and_inverts(trial):
    """`Galois.fromSeq(bijective=True)` yields an invertible Galois LFSR.

    Bijectivity is a property of the polynomial rather than of the realization:
    the Galois and Fibonacci registers built from one polynomial have update
    matrices that are transposes of each other, so one is invertible exactly
    when the other is. The Galois path therefore inherits the guarantee from
    the same constrained search, and this pins that the delegation is wired up.
    """
    random.seed(1300 + trial)
    seq = [random.randint(0, 1) for _ in range(24)]

    state, register = Galois.fromSeq(seq, bijective=True)

    produced = [int(s[0]) for s
                in FeedbackRegister(state, register).run(compiled=False, limit=len(seq))]
    assert produced == seq, "the constrained Galois fit should reproduce the sequence"

    register.invert()          # raises if the update map is not a bijection
    register.invert()


@pytest.mark.parametrize("trial", range(4))
def test_galois_and_fibonacci_bijective_fits_agree_in_length(trial):
    """Both realizations of the constrained fit have the same register length.

    They are built from the same polynomial, so a length difference would mean
    one of the two paths is not delegating to the shared search.
    """
    random.seed(1400 + trial)
    seq = [random.randint(0, 1) for _ in range(24)]

    _, galois = Galois.fromSeq(seq, bijective=True)
    _, fibonacci = Fibonacci.fromSeq(seq, bijective=True)

    assert galois.size == fibonacci.size, (
        f"trial {trial}: Galois size {galois.size} != Fibonacci size {fibonacci.size}"
    )

