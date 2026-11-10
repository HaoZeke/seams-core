#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Load the in-situ checks. A broken import fails this process."""

import importlib.util
import os
from pathlib import Path


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    root = Path(os.environ["SRC_SEAMS"])
    load(root / "tests" / "insitu" / "adapter_cov.py", "adapter_cov")
    load(root / "tests" / "insitu" / "test_shipped_in_dseams.py", "test_shipped")
    print("import checks exit 0", flush=True)


main()
