"""Progress reporting for long-running computations, built on `logging`.

Library code reports what it is doing through a `Reporter`, obtained with
`get_logger(__name__)`. Three calls describe the shape of a computation:

- `@log.stage("Offline phase")` on a function opens one level of nesting for
  the whole call and closes it -- with the call's total time -- when the
  function returns or raises.
- `log.step("Generating equations")` is a plain statement: it ends the
  current step of the enclosing stage (logging that step's time) and starts
  the next. The last step ends with the stage.
- `log.progress("Equations found", total=n)` returns a `Meter` to update from
  a loop. In a terminal it is drawn as a live bar; it ends as one summary line.

Everything else is ordinary logging (`log.info(...)`, `log.debug(...)`), and
every line is indented by how many stages and steps enclose it, so the output
reads as a tree of the computation::

    Offline phase (Reduced Algebraic Attack)
    |   Monomial profile
    |   Monomial profile finished: 0.42 s
    |   Generating equations
    |   |   Equations found: 52/52 in 0.37 s (140 eq/s)
    |   Generating equations finished: 0.39 s
    Offline phase (Reduced Algebraic Attack) finished: 0.83 s

Lines are printed once and never redrawn; only the bars of open meters are
live. Output goes to standard output, warnings and errors only by default;
`PyPR.logging.level = logging.INFO` shows a run's phases, results and progress, and
`"DEBUG"` everything. A stage, step or meter below the level costs almost
nothing and adds no nesting.
"""
import contextvars
import functools
import inspect
import logging
import time
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any, ParamSpec, TypeVar

ROOT_LOGGER = "PyPR"
INDENT = "|   "

# how often an open meter writes a progress line to non-live outputs (files, CI)
PROGRESS_LOG_INTERVAL = 10.0
# how often an open meter pushes its count to a live bar
_LIVE_PUSH_INTERVAL = 0.1

P = ParamSpec("P")
R = TypeVar("R")


# ── formatting helpers ───────────────────────────────────────────────────────

def format_duration(seconds: float) -> str:
    """Format a duration for a summary line: `0.42 s`, `12.3 s`, `3 min 12 s`, `1 h 03 min`.

    :param seconds: The duration.
    :type seconds: float
    :return: The formatted duration.
    :rtype: str
    """
    if seconds < 10:
        return f"{seconds:.2f} s"
    if seconds < 60:
        return f"{seconds:.1f} s"
    if seconds < 3600:
        minutes, secs = divmod(round(seconds), 60)
        return f"{minutes} min {secs:02d} s"
    hours, minutes = divmod(round(seconds) // 60, 60)
    return f"{hours} h {minutes:02d} min"


def _format_rate(count: float, seconds: float, unit: str) -> str:
    rate = count / seconds if seconds > 0 else 0.0
    for scale, suffix in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if rate >= scale:
            number = f"{rate / scale:.1f}{suffix}"
            break
    else:
        number = f"{rate:.0f}" if rate >= 10 else f"{rate:.2g}"
    return f"{number} {unit}/s" if unit else f"{number}/s"


def _format_count(count: int, total: int | None) -> str:
    return f"{count:,}/{total:,}" if total is not None else f"{count:,}"


# ── nesting state ────────────────────────────────────────────────────────────

@dataclass
class _Scope:
    """One call of a `stage`-decorated function, and the step open inside it."""
    logger: logging.Logger
    title: str
    level: int
    depth: int                      # indentation of the stage's own header
    start: float
    shown: bool                     # whether the stage's level is enabled
    step_title: str | None = None
    step_level: int = logging.INFO
    step_start: float = 0.0
    step_shown: bool = False
    meters: list["Meter"] = field(default_factory=list)

    @property
    def body_depth(self) -> int:
        """Indentation of lines logged in the stage outside any step."""
        return self.depth + (1 if self.shown else 0)

    @property
    def inner_depth(self) -> int:
        """Indentation of lines logged inside the current step."""
        return self.body_depth + (1 if self.step_title is not None and self.step_shown else 0)


# The open scopes, innermost last. A context variable, so nesting is tracked
# per thread (and per asyncio task) without being passed down explicitly.
_scopes: contextvars.ContextVar[tuple[_Scope, ...]] = contextvars.ContextVar(
    "pypr_scopes", default=()
)


def current_depth() -> int:
    """The indentation level a line logged now would get.

    :return: The number of enclosing, visible stages and steps.
    :rtype: int
    """
    scopes = _scopes.get()
    return scopes[-1].inner_depth if scopes else 0


# ── the reporter used by library code ────────────────────────────────────────

class Reporter:
    """A logger for one module, with stages, steps and meters.

    Obtain one with `get_logger(__name__)`. The plain logging methods
    (`debug`, `info`, `warning`, `error`) take the usual `%`-style arguments and
    indent the line to the current depth.

    :param name: The logger name, normally the module's `__name__`.
    :type name: str
    """

    def __init__(self, name: str):
        self.logger = logging.getLogger(name)

    # -- plain lines ---------------------------------------------------------
    def log(self, level: int, msg: str, *args: Any) -> None:
        """Log a line at `level`, indented to the current depth.

        :param level: The logging level.
        :type level: int
        :param msg: The message, with `%`-style placeholders.
        :type msg: str
        :param args: Values for the placeholders.
        :type args: Any
        """
        if self.logger.isEnabledFor(level):
            self.logger.log(level, msg, *args, extra={"pypr_depth": current_depth()}, stacklevel=3)

    def debug(self, msg: str, *args: Any) -> None:
        """Log a DEBUG line (see `log`)."""
        self.log(logging.DEBUG, msg, *args)

    def info(self, msg: str, *args: Any) -> None:
        """Log an INFO line (see `log`)."""
        self.log(logging.INFO, msg, *args)

    def warning(self, msg: str, *args: Any) -> None:
        """Log a WARNING line (see `log`)."""
        self.log(logging.WARNING, msg, *args)

    def error(self, msg: str, *args: Any) -> None:
        """Log an ERROR line (see `log`)."""
        self.log(logging.ERROR, msg, *args)

    def _emit(self, level: int, depth: int, msg: str, *args: Any, progress: bool = False) -> None:
        # lines whose depth is known already (headers, summaries): no lookup
        self.logger.log(
            level, msg, *args,
            extra={"pypr_depth": depth, "pypr_progress": progress},
            stacklevel=4,
        )

    # -- stages ----------------------------------------------------------------
    def stage(self,
        title: str,
        level: int = logging.INFO
    ) -> Callable[[Callable[P, R]], Callable[P, R]]:
        """Decorate a function so each call is one level of nesting in the output.

        The stage's title is logged when the call starts, everything logged
        during the call is indented one level beneath it, and a summary line
        with the call's total time is logged when it returns. If the call
        raises, the summary says so and the exception propagates unchanged.
        A stage below the enabled level logs nothing and adds no nesting, but
        its steps and lines are still logged at their own levels.

        :param title: The title of the stage.
        :type title: str
        :param level: The level of the stage's own lines. Defaults to INFO.
        :type level: int
        :raises TypeError: If applied to a generator or coroutine function,
            whose body does not run within the call.
        :return: The decorator.
        :rtype: Callable[[Callable[P, R]], Callable[P, R]]
        """
        def decorate(fn: Callable[P, R]) -> Callable[P, R]:
            if inspect.isgeneratorfunction(fn) or inspect.iscoroutinefunction(fn):
                raise TypeError(
                    f"stage cannot wrap {fn.__qualname__}: a generator or coroutine body "
                    f"runs after the call returns, outside the stage"
                )

            @functools.wraps(fn)
            def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
                scope = self._open_stage(title, level)
                token = _scopes.set(_scopes.get() + (scope,))
                try:
                    result = fn(*args, **kwargs)
                except BaseException as exc:
                    self._close_stage(scope, failure=exc)
                    raise
                finally:
                    _scopes.reset(token)
                self._close_stage(scope)
                return result

            return wrapper
        return decorate

    def _open_stage(self, title: str, level: int) -> _Scope:
        shown = self.logger.isEnabledFor(level)
        scope = _Scope(self.logger, title, level, current_depth(), time.perf_counter(), shown)
        if shown:
            self._emit(level, scope.depth, "%s", title)
        return scope

    def _close_stage(self, scope: _Scope, failure: BaseException | None = None) -> None:
        self._end_step(scope, failure)
        if not scope.shown:
            return
        elapsed = format_duration(time.perf_counter() - scope.start)
        if failure is None:
            self._emit(scope.level, scope.depth, "%s finished: %s", scope.title, elapsed)
        else:
            self._emit(max(scope.level, logging.WARNING), scope.depth, "%s failed after %s (%s)",
                       scope.title, elapsed, type(failure).__name__)

    # -- steps -----------------------------------------------------------------
    def step(self,
        title: str,
        level: int | None = None
    ) -> None:
        """End the current step of the enclosing stage, and start a new one.

        The step's title is logged now, lines logged until the next step are
        indented beneath it, and a summary line with its time is logged when
        it ends: at the next `step`, or when the stage returns. Meters opened
        during the step are closed with it. Outside any stage, a step is just
        its title line.

        :param title: The title of the step.
        :type title: str
        :param level: The level of the step's own lines. Defaults to the level
            of the enclosing stage (INFO outside any stage).
        :type level: int | None
        """
        scopes = _scopes.get()
        if not scopes:
            self.log(logging.INFO if level is None else level, "%s", title)
            return

        scope = scopes[-1]
        self._end_step(scope)
        scope.step_title = title
        scope.step_level = scope.level if level is None else level
        scope.step_start = time.perf_counter()
        scope.step_shown = scope.logger.isEnabledFor(scope.step_level)
        if scope.step_shown:
            self._emit(scope.step_level, scope.body_depth, "%s", title)

    def _end_step(self, scope: _Scope, failure: BaseException | None = None) -> None:
        for meter in list(scope.meters):
            meter.close(quiet=failure is not None)
        scope.meters.clear()

        if scope.step_title is None:
            return
        if scope.step_shown and failure is None:
            elapsed = format_duration(time.perf_counter() - scope.step_start)
            self._emit(scope.step_level, scope.body_depth, "%s finished: %s", scope.step_title, elapsed)
        scope.step_title = None

    # -- meters ----------------------------------------------------------------
    def progress(self,
        description: str,
        total: int | None = None,
        level: int | None = None,
        unit: str = ""
    ) -> "Meter":
        """Open a meter counting progress through a loop.

        In a live console the meter is a bar (or, without a total, a counter
        with a spinner) that updates in place; other outputs get a progress
        line at most every `PROGRESS_LOG_INTERVAL` seconds. When the meter
        closes -- explicitly, when its step or stage ends, or when its count
        reaches `total` -- it is replaced by one summary line such as
        `Equations found: 52/52 in 0.37 s (140 eq/s)`.

        :param description: What is being counted.
        :type description: str
        :param total: The count at completion, if known. Defaults to None.
        :type total: int | None
        :param level: The meter's level. Defaults to the level of the enclosing
            step or stage (INFO outside any stage).
        :type level: int | None
        :param unit: A unit for the rate, e.g. `"eq"` for `140 eq/s`. Defaults
            to none.
        :type unit: str
        :return: The meter.
        :rtype: Meter
        """
        scopes = _scopes.get()
        if level is None:
            if scopes and scopes[-1].step_title is not None:
                level = scopes[-1].step_level
            elif scopes:
                level = scopes[-1].level
            else:
                level = logging.INFO

        meter = Meter(self, description, total, level, unit)
        if scopes:
            scopes[-1].meters.append(meter)
        return meter


class Meter:
    """Progress through one loop. Create with `Reporter.progress`.

    Updating a meter whose level is not enabled only increments `count`, so
    meters can be updated freely from hot loops.

    :ivar count: The progress so far.
    :ivar total: The count at completion, if known.
    :ivar status: Free text shown after the count while the meter is open.
    """

    def __init__(self, reporter: Reporter, description: str, total: int | None, level: int, unit: str):
        self.reporter = reporter
        self.description = description
        self.total = total
        self.level = level
        self.unit = unit
        self.count = 0
        self.status = ""
        self.depth = current_depth()
        self.start = time.perf_counter()
        self.closed = False

        self.shown = reporter.logger.isEnabledFor(level)
        self._last_push = 0.0
        self._last_log = self.start
        self._displays: list[_LiveDisplay] = []
        if self.shown:
            for display in _live_displays_for(reporter.logger, level):
                display.add(self)
                self._displays.append(display)

    def update(self, amount: int = 1) -> None:
        """Advance the count by `amount`.

        :param amount: How much to add. Defaults to 1.
        :type amount: int
        """
        self.count += amount
        if self.shown:
            self._tick()

    def update_to(self, count: int) -> None:
        """Set the count directly, e.g. to a size that grows unevenly.

        :param count: The new count.
        :type count: int
        """
        self.count = count
        if self.shown:
            self._tick()

    def set_status(self, text: str) -> None:
        """Show `text` after the count while the meter is open.

        :param text: The status text.
        :type text: str
        """
        self.status = text
        if self.shown:
            self._tick()

    def _tick(self) -> None:
        if self.closed:
            return
        if self.total is not None and self.count >= self.total:
            self.close()
            return
        now = time.perf_counter()
        if self._displays and now - self._last_push >= _LIVE_PUSH_INTERVAL:
            self._last_push = now
            for display in self._displays:
                display.update(self)
        if now - self._last_log >= PROGRESS_LOG_INTERVAL:
            self._last_log = now
            rate = _format_rate(self.count, now - self.start, self.unit)
            status = f" -- {self.status}" if self.status else ""
            self.reporter._emit(self.level, self.depth, "%s: %s (%s)%s", self.description,
                                _format_count(self.count, self.total), rate, status, progress=True)

    def close(self,
        summary: str | None = None,
        quiet: bool = False
    ) -> None:
        """Close the meter: remove its bar, and log its summary line.

        Closing twice does nothing.

        :param summary: A line to log instead of the default
            `description: count/total in time (rate)`. Defaults to None.
        :type summary: str | None
        :param quiet: If True, log no summary line at all. Defaults to False.
        :type quiet: bool
        """
        if self.closed:
            return
        self.closed = True
        for display in self._displays:
            display.remove(self)
        if not self.shown or quiet:
            return
        if summary is None:
            elapsed = time.perf_counter() - self.start
            summary = (f"{self.description}: {_format_count(self.count, self.total)} "
                       f"in {format_duration(elapsed)} ({_format_rate(self.count, elapsed, self.unit)})")
        self.reporter._emit(self.level, self.depth, "%s", summary)


def get_logger(name: str) -> Reporter:
    """The `Reporter` for a module: `log = get_logger(__name__)`.

    :param name: The logger name, normally the module's `__name__`.
    :type name: str
    :return: A reporter wrapping `logging.getLogger(name)`.
    :rtype: Reporter
    """
    return Reporter(name)


# ── output ───────────────────────────────────────────────────────────────────

class _DepthFilter(logging.Filter):
    """Give every record the indentation string its formatter prints.

    Records from `Reporter` carry their depth; any other record logged inside
    a stage (for example through a plain `logging.getLogger`) is indented to
    the depth current when it is handled.
    """

    def filter(self, record: logging.LogRecord) -> bool:
        depth = getattr(record, "pypr_depth", None)
        if depth is None:
            depth = current_depth()
        record.pypr_indent = INDENT * depth
        return True


class _LiveDisplay:
    """The live bars at the bottom of a rich console, one per open meter.

    The bars exist only while some meter is open: the display starts with the
    first and stops (erasing them) with the last, so the console behaves
    normally between loops. Log lines printed through the same console appear
    above the bars.
    """

    def __init__(self, console: Any):
        self.console = console
        self.progress: Any = None
        self.tasks: dict[int, Any] = {}

    def add(self, meter: Meter) -> None:
        from rich.progress import (
            BarColumn,
            Progress,
            SpinnerColumn,
            TextColumn,
            TimeRemainingColumn,
        )

        if self.progress is None:
            self.progress = Progress(
                TextColumn("{task.fields[indent]}{task.description}", markup=False),
                SpinnerColumn(finished_text=""),
                BarColumn(bar_width=30),
                TextColumn("{task.fields[count]}", markup=False),
                TextColumn("{task.fields[rate]}", markup=False),
                TimeRemainingColumn(),
                TextColumn("{task.fields[status]}", markup=False),
                console=self.console,
                transient=True,
                refresh_per_second=10,
            )
            self.progress.start()
        self.tasks[id(meter)] = self.progress.add_task(
            meter.description, total=meter.total, completed=meter.count,
            indent=INDENT * meter.depth, count=_format_count(meter.count, meter.total),
            rate="", status=meter.status,
        )

    def update(self, meter: Meter) -> None:
        if id(meter) not in self.tasks:
            return
        elapsed = time.perf_counter() - meter.start
        self.progress.update(
            self.tasks[id(meter)], completed=meter.count,
            count=_format_count(meter.count, meter.total),
            rate=_format_rate(meter.count, elapsed, meter.unit), status=meter.status,
        )

    def remove(self, meter: Meter) -> None:
        if id(meter) not in self.tasks:
            return
        self.progress.remove_task(self.tasks.pop(id(meter)))
        if not self.tasks:
            self.progress.stop()
            self.progress = None


def _live_displays_for(logger: logging.Logger, level: int) -> list[_LiveDisplay]:
    """The live displays a meter on `logger` at `level` should draw on.

    Found the way logging itself routes a record -- up the logger hierarchy
    while records propagate -- so a console handler draws bars exactly when it
    would print the meter's lines, and detaching it with the standard
    `removeHandler` also stops its bars.
    """
    displays = []
    current: logging.Logger | None = logger
    while current is not None:
        for handler in current.handlers:
            if isinstance(handler, _ConsoleHandler) and handler.live is not None and level >= handler.level:
                displays.append(handler.live)
        if not current.propagate:
            break
        current = current.parent
    return displays


class _ConsoleHandler(logging.Handler):
    """Prints records through a rich console on standard output, above any live bars.

    The console is created the first time something is printed, so importing
    PyPR does not import rich. It writes to whatever `sys.stdout` is at each
    write, which is what makes `python run.py > out.txt` capture the output.
    Whether bars are drawn is decided when the console is created: on a
    terminal or in a notebook they are; on anything else (a pipe, a file) open
    meters write a progress line every `PROGRESS_LOG_INTERVAL` seconds instead.
    """

    def __init__(self, level: int):
        super().__init__(level)
        self._console: Any = None
        self._live: _LiveDisplay | None = None
        self.setFormatter(_ConsoleFormatter())
        self.addFilter(_DepthFilter())

    def use_console(self, console: Any) -> None:
        """Print through `console` (a `rich.console.Console`) from now on.

        :param console: The console, e.g. one with a fixed width or writing to
            a buffer in a test.
        :type console: rich.console.Console
        """
        self._console = console
        self._live = _LiveDisplay(console) if console.is_terminal or console.is_jupyter else None

    @property
    def console(self) -> Any:
        if self._console is None:
            from rich.console import Console
            self.use_console(Console())
        return self._console

    @property
    def live(self) -> "_LiveDisplay | None":
        if self._console is None:
            # creating the console is what decides whether bars are drawn
            _ = self.console
        return self._live

    def emit(self, record: logging.LogRecord) -> None:
        # a live console shows progress as bars; the periodic lines are for other outputs
        if self.live is not None and getattr(record, "pypr_progress", False):
            return
        try:
            self.console.print(self.format(record), markup=False, highlight=False, soft_wrap=True)
        # the logging.Handler contract: a handler never raises into the code
        # that logged; handleError reports the failure instead
        except Exception:  # noqa: BLE001
            self.handleError(record)


class _ConsoleFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        message = record.getMessage()
        if record.levelno >= logging.WARNING:
            message = f"{record.levelname}: {message}"
        return f"{getattr(record, 'pypr_indent', '')}{message}"


_LEVELS = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL,
    "OFF": logging.CRITICAL + 1,
}


class LoggingSettings:
    """PyPR's output settings, available as `PyPR.logging`::

        import logging
        import PyPR

        PyPR.logging.level = logging.INFO    # or the name: "INFO"

    The level decides how much prints, from least to most: `"OFF"`,
    `logging.ERROR`, `logging.WARNING` (the default), `logging.INFO` (phases,
    results and progress bars -- the normal view of a run) and
    `logging.DEBUG` (everything, including per-step solver and store detail).
    Output goes to standard output; pipe it to keep a record of a run.

    There is nothing else to set, and assigning any other attribute raises,
    so a misspelt setting fails loudly rather than being silently ignored.
    """

    __slots__ = ()

    @property
    def level(self) -> int:
        """The lowest level printed, as a `logging` level number.

        Set it to a `logging` level number or a level name (any case,
        including `"OFF"`).

        :raises ValueError: If set to an unknown level name.
        """
        return logging.getLogger(ROOT_LOGGER).level

    @level.setter
    def level(self, level: int | str) -> None:
        if isinstance(level, str):
            if level.upper() not in _LEVELS:
                raise ValueError(f"unknown log level {level!r}; use one of {', '.join(_LEVELS)}")
            level = _LEVELS[level.upper()]
        logging.getLogger(ROOT_LOGGER).setLevel(level)
        _handler.setLevel(level)

    def __repr__(self) -> str:
        names = {number: name for name, number in _LEVELS.items()}
        name = names.get(self.level, str(self.level))
        return f"<PyPR logging settings: level={name}>"


settings = LoggingSettings()


# PyPR's output: warnings and errors by default, more through PyPR.logging. The
# logger does not propagate, so a logging configuration in the calling script
# (basicConfig, say) does not print every line a second time.
_handler = _ConsoleHandler(logging.WARNING)
_root = logging.getLogger(ROOT_LOGGER)
_root.addHandler(_handler)
_root.setLevel(logging.WARNING)
_root.propagate = False
