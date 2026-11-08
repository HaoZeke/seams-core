#!/usr/bin/env python3
"""Drop checkout prefixes from compile text.

A prefix is any directory whose device and inode match the checkout
named on the command line. The checkout comes from git rev-parse.
"""

import os
import sys


def prefixes_for(root: str, text: str) -> list[str]:
    target = os.stat(root)
    found = set()
    flattened = text.replace("=", " ").replace(":", " ")
    for token in flattened.split():
        if not token.startswith("/"):
            continue
        path = token
        while path and path != "/":
            try:
                st = os.stat(path)
            except OSError:
                path = os.path.dirname(path)
                continue
            if (st.st_dev, st.st_ino) == (target.st_dev, target.st_ino):
                found.add(path)
                break
            path = os.path.dirname(path)
    return sorted(found, key=len, reverse=True)


def strip_text(root: str, text: str) -> str:
    for prefix in prefixes_for(root, text):
        text = text.replace(prefix + "/", "")
    return text


def main() -> None:
    root = sys.argv[1]
    sys.stdout.write(strip_text(root, sys.stdin.read()))


if __name__ == "__main__":
    main()
