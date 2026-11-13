#!/usr/bin/env python3
# Copyright (c) d-SEAMS developers
"""Execute tutorial code blocks from top to bottom.

The pages are the Org tutorials, the quick start, and the language
tutorials. Bash and Python blocks are programs. A text block or an
example block beside a bash program is output that program must
contain, one line at a time. C blocks are excerpts of
docs/examples/quickstart.c, which this test compiles when libseams.so
is next to the seams binary.

Bootstrap lines (meson, pixi, nix, conda, pip, spack, cargo, gfortran,
julia, Rscript, module) are part of the block. Pass --skip-bootstrap to
run only the analysis lines, against the seams binary you pass in.
That is what the Meson docs suite does, because the suite has already
built the binary.
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ORG_SRC = re.compile(
    r"#\+begin_src\s+(\S+)[^\n]*\n(.*?)#\+end_src",
    re.I | re.S,
)
ORG_EXAMPLE = re.compile(
    r"#\+begin_example\n(.*?)#\+end_example",
    re.I | re.S,
)
FENCE = re.compile(r"```(\w+)?\n(.*?)```", re.S)
BOOTSTRAP = re.compile(
    r"^(?:sudo\s+)?(?:meson|ninja|pixi|nix|conda|micromamba|pip|pip3|spack|"
    r"cargo|gfortran|julia|Rscript|module|source|lua|lua5\.4|luajit)\b"
)
SEAMS_LINE = re.compile(r"(?:^|&&|\n)\s*(?:(?:\./)?(?:bbdir|build|result/bin)/)?seams\b")


def pages():
    found = sorted((ROOT / "docs" / "orgmode" / "tutorials").glob("*.org"))
    found = [path for path in found if path.name != "index.org"]
    found.append(ROOT / "docs" / "orgmode" / "quickstart.org")
    found.append(ROOT / "docs" / "quickstart.md")
    found.extend(sorted((ROOT / "docs" / "tutorials").glob("*.md")))
    return found


def org_blocks(text):
    marks = []
    for match in ORG_SRC.finditer(text):
        marks.append((match.start(), "code", match.group(1).lower(), match.group(2).strip("\n")))
    for match in ORG_EXAMPLE.finditer(text):
        marks.append((match.start(), "expect", "text", match.group(1).strip("\n")))
    marks.sort()
    return [(kind, lang, body) for _pos, kind, lang, body in marks]


def md_blocks(text):
    found = []
    for match in FENCE.finditer(text):
        lang = (match.group(1) or "").lower()
        body = match.group(2).strip("\n")
        if lang in {"text", "output"}:
            found.append(("expect", "text", body))
        elif lang:
            found.append(("code", lang, body))
    return found


def blocks_of(path):
    text = path.read_text(encoding="utf-8")
    if path.suffix == ".org":
        return org_blocks(text)
    return md_blocks(text)


def is_bootstrap(line):
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return False
    return BOOTSTRAP.match(stripped) is not None


def analysis_script(body):
    """Drop bootstrap lines. Keep seams commands and shell glue."""
    kept = []
    for line in body.splitlines():
        if is_bootstrap(line):
            continue
        kept.append(line)
    script = "\n".join(kept).strip()
    if not script or not SEAMS_LINE.search(script):
        return ""
    return script


def rewrite_seams(script, binary):
    """Replace the seams command, not a filename such as seams.env."""

    def swap(match):
        return match.group(1) + match.group(2) + str(binary)

    return re.sub(
        r"(^|&&|\|\||[;\n])(\s*)(?:\./(?:bbdir|build|result/bin)/)?seams(?=\s|$)",
        swap,
        script,
        flags=re.M,
    )


def find_binary(argv_path):
    if argv_path:
        path = Path(argv_path)
        if path.is_file():
            return path.resolve()
        raise SystemExit("seams binary not found: {}".format(argv_path))
    env = os.environ.get("SEAMS", "")
    candidates = []
    if env:
        candidates.append(Path(env))
    candidates.extend(
        [
            ROOT / "bbdir" / "seams",
            ROOT / "bbdir-eessi" / "seams",
            ROOT / "build" / "seams",
        ]
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise SystemExit("pass the seams binary, or set SEAMS")


def run_shell(script, binary, timeout):
    env = os.environ.copy()
    lib_dirs = [str(binary.parent), str(binary.parent / "src"), str(binary.parent / "bindings")]
    env["PATH"] = str(binary.parent) + os.pathsep + env.get("PATH", "")
    env["LD_LIBRARY_PATH"] = os.pathsep.join(lib_dirs + [env.get("LD_LIBRARY_PATH", "")])
    env["OMP_NUM_THREADS"] = env.get("OMP_NUM_THREADS", "1")
    done = subprocess.run(
        ["bash", "-e", "-c", rewrite_seams(script, binary)],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )
    return done.returncode, done.stdout, done.stderr


def run_python(body, namespace, timeout):
    """Execute one Python block in the page namespace. timeout is unused; exec is in-process."""
    del timeout
    import io
    from contextlib import redirect_stderr, redirect_stdout

    out = io.StringIO()
    err = io.StringIO()
    try:
        with redirect_stdout(out), redirect_stderr(err):
            exec(compile(body, "<tutorial>", "exec"), namespace, namespace)
    except Exception as exc:
        return 1, out.getvalue(), "{}\n{}".format(err.getvalue(), exc)
    return 0, out.getvalue(), err.getvalue()


def expect_lines(body):
    lines = []
    for line in body.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if "(" in stripped:
            return []
        lines.append(stripped)
    return lines


def check_expect(stdout, expected):
    missing = [line for line in expected if line not in stdout]
    return missing


def compile_quickstart(binary):
    lib = binary.parent / "src" / "libseams.so"
    header = ROOT / "include" / "seams.h"
    source = ROOT / "docs" / "examples" / "quickstart.c"
    if not lib.is_file() or not header.is_file() or not source.is_file():
        return None
    out = binary.parent / "docs_quickstart_c"
    cmd = [
        "cc",
        "-I{}".format(ROOT / "include"),
        str(source),
        str(lib),
        "-lm",
        "-o",
        str(out),
    ]
    done = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=False)
    if done.returncode != 0:
        return "compile failed: {}".format(done.stderr.strip()[:400])
    run = subprocess.run(
        [str(out)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
        env={**os.environ, "LD_LIBRARY_PATH": str(lib.parent) + os.pathsep + os.environ.get("LD_LIBRARY_PATH", "")},
    )
    if run.returncode != 0 or "cubic 8" not in run.stdout:
        return "quickstart.c printed {!r}".format(run.stdout.strip())
    return None


def main(argv):
    parser = argparse.ArgumentParser(description="Run tutorial code blocks")
    parser.add_argument("seams", nargs="?", help="Path to the seams binary")
    parser.add_argument(
        "--skip-bootstrap",
        action="store_true",
        help="Run analysis lines only. The binary is already built.",
    )
    parser.add_argument("--timeout", type=int, default=180)
    args = parser.parse_args(argv)
    binary = find_binary(args.seams)
    failures = []
    ran = 0
    for path in pages():
        namespace = {"__name__": "tutorial"}
        pending_expect = []
        last = None
        for kind, lang, body in blocks_of(path):
            rel = path.relative_to(ROOT)
            if kind == "expect":
                expected = expect_lines(body)
                if last and last[0] in {"bash", "sh", "shell", "console"} and expected:
                    missing = check_expect(last[1], expected)
                    if missing:
                        failures.append("{} missing {}\n{}".format(last[2], missing, last[1][:500]))
                else:
                    pending_expect = expected
                continue
            if lang in {"c", "cpp", "cxx", "lua"}:
                continue
            if lang not in {"bash", "sh", "shell", "console", "python", "py"}:
                continue
            script = body
            if lang in {"bash", "sh", "shell", "console"}:
                if args.skip_bootstrap:
                    script = analysis_script(body)
                    if not script:
                        pending_expect = []
                        last = None
                        continue
                code, stdout, stderr = run_shell(script, binary, args.timeout)
            else:
                first = body.strip().splitlines()[0] if body.strip() else ""
                if args.skip_bootstrap and is_bootstrap(first):
                    pending_expect = []
                    last = None
                    continue
                code, stdout, stderr = run_python(body, namespace, args.timeout)
            ran += 1
            label = "{} [{}]".format(rel, lang)
            if code != 0:
                failures.append("{}\n{}".format(label, (stderr or stdout)[:500]))
                pending_expect = []
                last = None
                continue
            if pending_expect and lang in {"bash", "sh", "shell", "console"}:
                missing = check_expect(stdout, pending_expect)
                if missing:
                    failures.append("{} missing {}\n{}".format(label, missing, stdout[:500]))
            pending_expect = []
            last = (lang, stdout, label)
    quick = compile_quickstart(binary)
    if quick:
        failures.append(quick)
    if failures:
        print("\n".join(failures), file=sys.stderr)
        print("{} tutorial blocks failed".format(len(failures)), file=sys.stderr)
        return 1
    print("tutorial blocks ok", ran)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
