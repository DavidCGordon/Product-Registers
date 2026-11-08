"""Sphinx configuration for the PyPR documentation.

Build from the repository root with::

    pip install -e ".[docs]"
    sphinx-build -b html docs docs/_build/html

The prose in ``theory/``, ``architecture/`` and ``conventions/`` is Markdown,
read through MyST. The API reference under ``api/`` is generated from the
docstrings on every build and is not checked in.

Hosted docs are built only from release and candidate tags. Each build's title
names its tag, and links that leave ``docs/`` point at the commit being built,
so the docs for a release never send readers to newer code. See "Documentation
versions" in ``Branching and Releases.md``.
"""
import os
import re
import subprocess
import tomllib
from pathlib import Path
from urllib.parse import quote, unquote

DOCS = Path(__file__).resolve().parent


def built_ref():
    """Name the git branch or tag this build is made from.

    Read the Docs states it in ``READTHEDOCS_GIT_IDENTIFIER``; a local build
    asks git for the current branch, and falls back to ``main`` outside a clone.

    :return: the branch or tag name
    :rtype: str
    """
    if "READTHEDOCS_GIT_IDENTIFIER" in os.environ:
        return os.environ["READTHEDOCS_GIT_IDENTIFIER"]
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            cwd=DOCS, capture_output=True, text=True, check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return "main"
    return result.stdout.strip()


REF = built_ref()
# A pull-request preview on Read the Docs is named by its PR number, which is no
# git ref; its commit hash is, and is exact for every other build too.
REPO_URL = (
    "https://github.com/DavidCGordon/Product-Registers/blob/"
    + os.environ.get("READTHEDOCS_GIT_COMMIT_HASH", REF) + "/"
)

project = "PyPR"
author = "David Gordon"
copyright = "David Gordon"
# pyproject.toml is the only place the version is written; Sphinx shows it in
# the API reference's metadata.
release = version = tomllib.loads((DOCS.parent / "pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]

extensions = [
    "myst_parser",
    "sphinx.ext.apidoc",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
]

# The documents write math as $...$ and $$...$$, and GitHub renders that too.
myst_enable_extensions = ["dollarmath", "amsmath"]
# Headings down to this level get anchors, so links of the form Doc.md#section resolve.
myst_heading_anchors = 4

exclude_patterns = ["_build"]

# One page per module, generated from src/. Several packages have no
# __init__.py (Cryptanalysis/Attacks, Tools/AlgClosure), hence the namespaces.
apidoc_modules = [
    {
        "path": "../src/PyPR",
        "destination": "api",
        "separate_modules": True,
        "module_first": True,
        "implicit_namespaces": True,
        # A package that re-exports a class would otherwise document it a
        # second time; each object appears on the page of the module defining it.
        "automodule_options": {"members", "undoc-members", "show-inheritance", "ignore-module-all"},
    },
]
autodoc_member_order = "bysource"
# The docstrings mark code with single backticks, as Markdown does.
default_role = "code"
autodoc_typehints = "description"



def link_outside_docs(app, docname, source):
    """Point links that leave ``docs/`` at the repository.

    The documents link to tests, agent files and the README by relative path,
    which works on GitHub but names nothing inside the built site. Images are
    left alone: Sphinx copies those in from wherever they are.
    """
    folder = (DOCS / docname).parent

    def link(match):
        target = (folder / unquote(match.group(2))).resolve()
        if DOCS in target.parents or not target.exists():
            return match.group(0)
        return match.group(1) + "(" + REPO_URL + quote(target.relative_to(DOCS.parent).as_posix()) + (match.group(3) or "") + ")"

    source[0] = re.sub(r"((?<!!)\[[^\]\n]*\])\((\.\.[^()#\s]*)(#[^()\s]*)?\)", link, source[0])


def unwrap_environments(app, docname, source):
    """Take the ``$$`` off display math that is a single AMS environment.

    GitHub needs ``$$`` around ``\\begin{align*}``. Sphinx wraps any ``$$``
    block containing a line break in a ``split`` environment, and ``align``
    inside ``split`` is an error in MathJax ("Erroneous nesting of equation
    structures"). Without the ``$$``, MyST's ``amsmath`` extension passes the
    environment through as it stands.
    """
    source[0] = re.sub(
        r"^([ \t]*)\$\$[ \t]*\n\s*(\\begin\{(align|gather|multline|alignat|flalign|equation)\*?\}.*?\\end\{\3\*?\})\s*\n[ \t]*\$\$[ \t]*$",
        r"\1\2",
        source[0],
        flags=re.DOTALL | re.MULTILINE,
    )


def skip_reexports(app, what, name, obj, skip, options):
    """Document each object once, on the page of the module that defines it.

    ``PyPR`` re-exports the ``FeedbackRegister`` class, whose name is also that
    of the module it lives in, so documenting it under the package as well
    gives two objects the address ``PyPR.FeedbackRegister``.
    """
    current = app.env.temp_data.get("autodoc:module")
    if what == "module" and current and getattr(obj, "__module__", current) != current:
        return True
    return skip


def setup(app):
    app.connect("source-read", link_outside_docs)
    app.connect("source-read", unwrap_environments)
    app.connect("autodoc-skip-member", skip_reexports)


html_theme = "furo"
# Hosted docs are built only from version tags: releases (v1.4.0) and
# candidates (v1.4.0-rc2), each named by its tag. Branch names never start with
# "v". main always holds exactly the newest release, so it is named by the
# version in its pyproject.toml. Any other build, a local one on dev for
# instance, describes unreleased work, whatever version it last recorded.
if REF.startswith("v"):
    html_title = f"PyPR {REF[1:]}"
elif REF == "main":
    html_title = f"PyPR {release}"
else:
    html_title = "PyPR (unreleased)"
