#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Load the in-situ checks and show each negative helper can fail."""

import importlib.util
import os
from pathlib import Path


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def show_fail(fn, needle):
    try:
        fn()
    except SystemExit as err:
        text = err.code if isinstance(err.code, str) else str(err)
        if needle not in text:
            raise SystemExit("raised %r" % text)
        print("fail-shown", needle, flush=True)
        return
    raise SystemExit("did not fail: " + needle)


def main():
    root = Path(os.environ["SRC_SEAMS"])
    adapter = load(root / "tests" / "insitu" / "adapter_cov.py", "adapter_cov")
    shipped = load(root / "tests" / "insitu" / "test_shipped_in_dseams.py", "test_shipped")
    show_fail(lambda: adapter.must_raise(lambda: None, "set "), "did not raise")
    show_fail(
        lambda: adapter.must_raise(lambda: (_ for _ in ()).throw(SystemExit("other")), "set "),
        "raised",
    )
    show_fail(lambda: shipped.assert_reject(lambda: None, "set X"), "reject did not fire")
    show_fail(
        lambda: shipped.assert_reject(lambda: (_ for _ in ()).throw(SystemExit("nope")), "other"),
        "reject was",
    )
    show_fail(lambda: shipped.need("DSEAMS_NOT_SET"), "set DSEAMS_NOT_SET")

    def caught(fn, needle):
        try:
            show_fail(fn, needle)
        except SystemExit as exc:
            return str(exc)
        return ""

    got_ok = caught(lambda: shipped.need("DSEAMS_NOT_SET"), "set ")
    got_miss = caught(lambda: shipped.need("PATH"), "nope")
    got_needle = caught(lambda: shipped.need("DSEAMS_NOT_SET"), "not-the-needle")
    def require_arms(arms, want):
        if arms != want:
            raise SystemExit("show_fail arms %s" % arms)

    arms = [got_ok == "", got_miss.startswith("did not fail"), got_needle.startswith("raised ")]
    def keep(exc, needle):
        if needle not in str(exc):
            raise SystemExit("bad arm")

    require_arms(arms, [True, True, True])
    try:
        require_arms([False], [True])
    except SystemExit as exc:
        keep(exc, "show_fail arms")
        try:
            keep(exc, "absent-needle")
        except SystemExit:
            pass
    print("import checks exit 0", flush=True)


main()
