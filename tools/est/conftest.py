"""Planted-bug support for the estimator tests.

mutate(module, name, old, new) recompiles a copy of module.name with the source
text `old` replaced by `new` (which must occur exactly once, so a refactor that
removes the mutated line fails loudly instead of silently testing nothing) and
monkeypatches it in for the duration of one test. The copy runs in the
module's globals, so calls into other module functions see any other patches.
"""
from __future__ import annotations
import __future__

import inspect
import textwrap

import pytest


@pytest.fixture
def mutate(monkeypatch):
    def apply(module, name: str, old: str, new: str) -> None:
        src = textwrap.dedent(inspect.getsource(getattr(module, name)))
        n = src.count(old)
        if n != 1:
            raise AssertionError(f"mutation target occurs {n} times in {module.__name__}.{name}: {old!r}")
        code = compile(src.replace(old, new), f"<mutant {module.__name__}.{name}>", "exec",
                       flags=__future__.annotations.compiler_flag, dont_inherit=True)
        ns: dict = {}
        exec(code, module.__dict__, ns)
        monkeypatch.setattr(module, name, ns[name])

    return apply
