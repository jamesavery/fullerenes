"""Regression tests for BuckyGen child-process teardown via the binding.

These guard the long-lived-host scenario the teardown rework targets: a Python
interpreter that opens/closes many generators and may have installed its own
signal handlers. The headline guard is that close() must not hang when a host
SIGTERM handler is installed (the forked child resets SIGTERM to SIG_DFL).

Each potentially-hanging call is bounded by SIGALRM so a regression FAILS the
test instead of hanging the suite.
"""
import os
import shutil
import signal
import subprocess

import pytest

import fullerenes as fl


class _Timeout(Exception):
    pass


def _bounded(seconds, fn):
    """Run fn() with a SIGALRM deadline (pytest runs in the main thread)."""
    def _raise(*_):
        raise _Timeout()
    old = signal.signal(signal.SIGALRM, _raise)
    signal.alarm(seconds)
    try:
        return fn()
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def _zombie_children():
    if shutil.which("ps") is None:
        pytest.skip("ps not available")
    out = subprocess.run(["ps", "-o", "stat=", "--ppid", str(os.getpid())],
                         capture_output=True, text=True).stdout
    return [s for s in out.split() if s.startswith("Z")]


def test_close_does_not_hang_with_host_sigterm_handler():
    # A host SIGTERM handler must not make close() hang: the child resets SIGTERM
    # to SIG_DFL, so killpg terminates it and waitpid returns. Without the reset,
    # killpg would run the inherited (non-terminating) handler in the child and
    # close()'s waitpid would block forever -> _Timeout here.
    old = signal.signal(signal.SIGTERM, lambda *_: None)
    try:
        def body():
            for _ in range(10):
                g = fl.buckygen(60)
                next(iter(g))
                g.close()
        _bounded(30, body)
    finally:
        signal.signal(signal.SIGTERM, old)


def test_closing_one_generator_does_not_kill_siblings():
    # Per-queue process group: closing one generator must not SIGTERM another's
    # child (the old killpg-whole-group behaviour did).
    a = fl.buckygen(60)
    b = fl.buckygen(70)
    next(iter(a))
    next(iter(b))
    a.close()
    assert _bounded(15, lambda: next(iter(b))).N == 37   # C70 dual: 37 vertices
    b.close()


def test_repeated_generators_leave_no_zombies():
    # waitpid in stop() must reap each child -> no zombie accumulation.
    for i in range(25):
        g = fl.buckygen(60 + 2 * (i % 6))
        next(iter(g))
        g.close()
    assert _zombie_children() == []


def test_double_close_is_a_safe_noop():
    g = fl.buckygen(60)
    next(iter(g))
    g.close()
    g.close()           # idempotent: must not raise or kill a reused pid
    assert _zombie_children() == []
