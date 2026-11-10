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
    show_fail(lambda: show_fail(lambda: None, "nope"), "did not fail")
    show_fail(
        lambda: show_fail(lambda: (_ for _ in ()).throw(SystemExit("other")), "needle"),
        "raised",
    )
    print("import checks exit 0", flush=True)


main()
