"""Tests for PyPR.Reporting: nested stages, steps and progress meters over logging.

Every line is checked as text, because the text is the product: indentation is
how the output shows which work happened inside which, and the summary lines
are where the timings and counts end up. Timings vary, so they are replaced by
a placeholder before comparing.
"""
import contextlib
import io
import logging
import re
import threading

import pytest
from rich.console import Console

import PyPR
from PyPR import Reporting
from PyPR.Reporting import format_duration, get_logger

log = get_logger("PyPR.test_reporting")


@contextlib.contextmanager
def captured(level: int | str = "INFO", terminal: bool = False):
    """Point PyPR's output at a buffer at `level`; yield a function that returns
    what was printed, timings masked, as a list of lines. A terminal console
    draws live bars; any other (like a pipe) gets periodic progress lines."""
    buffer = io.StringIO()
    handler = Reporting._handler
    previous = (logging.getLogger("PyPR").level, handler._console, handler._live)
    handler.use_console(Console(file=buffer, force_terminal=terminal, width=120))
    PyPR.logging.level = level

    def lines():
        # strip terminal control codes (the live bars) and mask timings/rates
        text = re.sub(r"\x1b\[[0-9;?]*[A-Za-z]", "", buffer.getvalue())
        text = re.sub(r"\d+(\.\d+)? (s|min|h)\b", "T", text)
        text = re.sub(r"\([\d.]+[kMG]? ?[a-z]*/s\)", "(R)", text)
        return [line for line in text.replace("\r", "\n").splitlines() if line.strip()]

    try:
        yield lines
    finally:
        PyPR.logging.level = previous[0]
        handler._console, handler._live = previous[1], previous[2]


# ── stages and steps ─────────────────────────────────────────────────────────

@log.stage("Outer")
def outer():
    log.step("First step")
    log.info("inside the first step")
    log.step("Second step")
    inner()
    log.info("after the inner stage")

@log.stage("Inner")
def inner():
    log.info("inside the inner stage")

def test_stages_and_steps_nest_and_report_their_times():
    with captured() as lines:
        outer()
    assert lines() == [
        "Outer",
        "|   First step",
        "|   |   inside the first step",
        "|   First step finished: T",
        "|   Second step",
        "|   |   Inner",
        "|   |   |   inside the inner stage",
        "|   |   Inner finished: T",
        "|   |   after the inner stage",
        "|   Second step finished: T",
        "Outer finished: T",
    ]

def test_a_stage_below_the_enabled_level_adds_no_nesting():
    """Its own lines are hidden, but lines it logs at an enabled level still
    appear -- at the depth of whatever encloses the hidden stage."""
    @log.stage("Detail", level=logging.DEBUG)
    def detail():
        log.step("Detail step")
        log.info("an INFO line inside a DEBUG stage")

    @log.stage("Visible")
    def visible():
        detail()

    with captured("INFO") as lines:
        visible()
    assert lines() == [
        "Visible",
        "|   an INFO line inside a DEBUG stage",
        "Visible finished: T",
    ]

    with captured("DEBUG") as lines:
        visible()
    assert lines()[1:4] == ["|   Detail", "|   |   Detail step", "|   |   |   an INFO line inside a DEBUG stage"]

def test_a_failing_stage_reports_it_and_the_exception_propagates():
    @log.stage("Fragile")
    def fragile():
        log.step("About to fail")
        meter = log.progress("Items", total=10)
        meter.update(3)
        raise ValueError("boom")

    with captured() as lines, pytest.raises(ValueError, match="boom"):
        fragile()
    # no step summary and no meter summary: neither finished
    assert lines() == ["Fragile", "|   About to fail", "WARNING: Fragile failed after T (ValueError)"]

def test_stage_refuses_a_generator_function():
    with pytest.raises(TypeError, match="generator or coroutine"):
        @log.stage("Lazy")
        def lazy():
            yield 1

def test_a_step_outside_any_stage_is_just_a_line():
    with captured() as lines:
        log.step("Loose step")
        log.info("next")
    assert lines() == ["Loose step", "next"]

def test_plain_loggers_inside_a_stage_are_indented_too():
    """A module logging through logging.getLogger rather than a Reporter still
    lands at the right depth, because the handler computes it."""
    plain = logging.getLogger("PyPR.test_reporting.plain")

    @log.stage("Wrapper")
    def wrapper():
        plain.info("from a plain logger")

    with captured() as lines:
        wrapper()
    assert lines()[1] == "|   from a plain logger"

def test_nesting_is_tracked_per_thread():
    """A stage open in one thread does not indent lines logged in another."""
    opened, release = threading.Event(), threading.Event()

    @log.stage("Background")
    def background():
        opened.set()
        release.wait(5)

    with captured() as lines:
        thread = threading.Thread(target=background)
        thread.start()
        opened.wait(5)
        log.info("main thread, while the background stage is open")
        release.set()
        thread.join()
    assert "main thread, while the background stage is open" in lines()


# ── meters ───────────────────────────────────────────────────────────────────

def test_meter_summary_keeps_the_count_out_of_the_total():
    @log.stage("Loop")
    def loop():
        found = log.progress("Equations found", total=52, unit="eq")
        for _ in range(37):
            found.update()
        found.close()
        log.info("done")

    with captured() as lines:
        loop()
    assert lines() == ["Loop", "|   Equations found: 37/52 in T (R)", "|   done", "Loop finished: T"]

def test_meter_closes_itself_at_its_total_and_ignores_later_updates():
    meter = None

    @log.stage("Loop")
    def loop():
        nonlocal meter
        meter = log.progress("Items", total=5)
        for _ in range(7):
            meter.update()
        log.info("after")

    with captured() as lines:
        loop()
    assert lines()[1:3] == ["|   Items: 5/5 in T (R)", "|   after"]
    assert meter is not None
    assert meter.closed

def test_meter_without_a_total_counts_and_takes_a_custom_summary():
    with captured() as lines:
        branch = log.progress("x_17 = 0")
        branch.update_to(120)
        branch.set_status("Basis: 30")
        branch.close(summary="x_17 = 0) found inconsistent")
        branch.close(summary="closing twice does nothing")
        quiet = log.progress("x_17 = 1")
        quiet.close(quiet=True)
    assert lines() == ["x_17 = 0) found inconsistent"]

def test_open_meters_close_with_their_step():
    @log.stage("Stage")
    def stage():
        log.step("Counting")
        meter = log.progress("Items")
        meter.update(4)
        log.step("Next")

    with captured() as lines:
        stage()
    assert lines()[1:4] == ["|   Counting", "|   |   Items: 4 in T (R)", "|   Counting finished: T"]

def test_meters_below_the_enabled_level_only_count():
    with captured("INFO") as lines:
        meter = log.progress("Hidden", total=3, level=logging.DEBUG)
        meter.update(2)
        meter.close()
    assert meter.count == 2
    assert not meter.shown
    assert lines() == []

def test_live_console_draws_bars_and_still_prints_the_summaries():
    """A terminal console runs the rich live display; the permanent lines are
    the same as for any other output."""
    @log.stage("Split")
    def split():
        branch_0, branch_1 = log.progress("x_3 = 0"), log.progress("x_3 = 1")
        branch_0.set_status("Processed: 10")
        branch_1.set_status("Processed: 12")
        branch_0.close(quiet=True)
        branch_1.close(quiet=True)
        guesses = log.progress("Guesses", total=8)
        for _ in range(8):
            guesses.update()

    with captured(terminal=True) as lines:
        split()
    permanent = lines()
    assert "Split" in permanent
    assert "|   Guesses: 8/8 in T (R)" in permanent
    assert "Split finished: T" in permanent

@pytest.mark.parametrize("terminal", [False, True], ids=["piped", "terminal"])
def test_piped_output_gets_periodic_progress_lines_instead_of_bars(terminal, monkeypatch):
    """Piped output cannot redraw a bar, so an open meter writes a line every
    PROGRESS_LOG_INTERVAL seconds there (here: every update). A terminal draws
    the bar and gets no such lines; both end with the same summary."""
    monkeypatch.setattr(Reporting, "PROGRESS_LOG_INTERVAL", 0.0)
    with captured(terminal=terminal) as lines:
        meter = log.progress("Items", total=100)
        for _ in range(3):
            meter.update()
        meter.close()
    progress_lines = [line for line in lines() if line.startswith("Items: ") and "/100 (" in line]
    assert len(progress_lines) == (0 if terminal else 3)
    assert lines()[-1] == "Items: 3/100 in T (R)"


# ── setup ────────────────────────────────────────────────────────────────────

def test_by_default_warnings_show_and_info_does_not():
    with captured(level=logging.getLogger("PyPR").level) as lines:
        outer()
        log.warning("something worth knowing")
    assert lines() == ["WARNING: something worth knowing"]

def test_level_off_silences_warnings_too():
    with captured(level="OFF") as lines:
        log.warning("not shown")
        log.error("not shown either")
    assert lines() == []

def test_level_rejects_unknown_names():
    with pytest.raises(ValueError, match="unknown log level 'LOUD'"):
        PyPR.logging.level = "LOUD"

def test_level_reads_back_as_a_logging_number():
    previous = PyPR.logging.level
    try:
        PyPR.logging.level = logging.DEBUG
        assert PyPR.logging.level == logging.DEBUG
        PyPR.logging.level = "info"
        assert PyPR.logging.level == logging.INFO
        assert repr(PyPR.logging) == "<PyPR logging settings: level=INFO>"
    finally:
        PyPR.logging.level = previous

def test_a_misspelt_setting_raises():
    with pytest.raises(AttributeError):
        PyPR.logging.levle = logging.INFO  # type: ignore[attr-defined]

def test_a_logging_config_in_the_calling_script_does_not_duplicate_lines(capsys):
    """PyPR's logger does not propagate, so a root handler set up by the caller
    (basicConfig, say) does not print PyPR's lines a second time."""
    root = logging.getLogger()
    extra = logging.StreamHandler()
    root.addHandler(extra)
    try:
        with captured() as lines:
            log.info("printed once")
    finally:
        root.removeHandler(extra)
    assert lines() == ["printed once"]
    assert "printed once" not in capsys.readouterr().err

@pytest.mark.parametrize(("seconds", "text"), [
    (0.4214, "0.42 s"), (12.34, "12.3 s"), (192, "3 min 12 s"), (3780, "1 h 03 min"),
])
def test_format_duration(seconds, text):
    assert format_duration(seconds) == text
