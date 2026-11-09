"""The coverage scorer rejects a directory, a mismatched source, and a missing path."""

import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def load_scorer():
    path = Path(__file__).resolve().parents[2] / "scripts" / "insitu_coverage.py"
    spec = importlib.util.spec_from_file_location("insitu_coverage", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run(module, root, *inputs):
    sys.argv = ["insitu_coverage.py", str(root), *[str(item) for item in inputs]]
    held = sys.stdout
    buffer = []

    class Grab:
        def write(self, text):
            buffer.append(text)

        def flush(self):
            return None

    grab = Grab()
    grab.flush()
    sys.stdout = grab
    try:
        try:
            module.main()
            code = 0
        except SystemExit as exc:
            code = exc.code
    finally:
        sys.stdout = held
    return code, "".join(buffer)


def need(cond, message):
    if not cond:
        raise SystemExit(message)


def main():
    module = load_scorer()
    try:
        need(False, "probe")
    except SystemExit:
        pass
    need(1 in module.physical("keep /* c */ tail\n"), "c block kept the tail")
    need(3 in module.physical("/* open\nstill\n*/ after\n"), "c block kept the line after it")
    need(module.physical("/* never ends\n") == set(), "an open c comment kept a line")
    need(1 in module.physical("head /* cut\n"), "c comment ate the code before it")
    need(1 in module.physical("code // note\n"), "c line comment ate the code")
    py = 'note = "/* never ends\\n"\nvalue = 1\n'
    need(module.physical(py, ".py") == {1, 2}, "a python string opened a c comment")
    quoted = '"""/* never\nends"""\nvalue = 1\n'
    need(module.physical(quoted, ".py") == {1, 2, 3}, "a python triple quote opened a c comment")
    with tempfile.TemporaryDirectory() as raw:
        root = Path(raw)
        git_env = os.environ.copy()
        git_env.update({
            "GIT_AUTHOR_NAME": "Rohit Goswami",
            "GIT_AUTHOR_EMAIL": "rgoswami@ieee.org",
            "GIT_COMMITTER_NAME": "Rohit Goswami",
            "GIT_COMMITTER_EMAIL": "rgoswami@ieee.org",
        })
        subprocess.run(["git", "init"], cwd=root, check=True, capture_output=True)
        (root / "anchor.txt").write_text("anchor\n")
        subprocess.run(["git", "add", "anchor.txt"], cwd=root, check=True)
        subprocess.run(
            ["git", "-c", "core.hooksPath=/dev/null", "commit", "-m", "anchor"],
            cwd=root, check=True, capture_output=True, env=git_env,
        )
        anchor = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
        saved_base = os.environ.get("DSEAMS_COVERAGE_BASE")
        os.environ["DSEAMS_COVERAGE_BASE"] = anchor
        folder = root / "notes"
        folder.mkdir()
        code, text = run(module, root, folder)
        need(code == 1 and "0 0 0.000000" in text, "directory row was not a failure")

        source = root / "widget.cpp"
        source.write_text("int other;\n")
        gcov = root / "widget.cpp.gcov"
        gcov.write_text("        -:    0:Source:widget.cpp\n        1:    1:int real;\n")
        code, text = run(module, root, gcov)
        need(code == 1 and "widget.cpp 0 0 0.000000" in text, "mismatched gcov was accepted")

        nested = root / "src"
        nested.mkdir()
        (nested / "widget.cpp").write_text("int real;\n")
        other = root / "other"
        other.mkdir()
        (other / "widget.cpp").write_text("int other;\n")
        code, text = run(module, root, gcov)
        need(code == 0 and "src/widget.cpp 1 1 1.000000" in text, "matching gcov was rejected")

        payload = {"files": {str(root / "absent.py"): {"executed_lines": [1], "missing_lines": []}}}
        js = root / "cov.json"
        js.write_text(json.dumps(payload))
        code, text = run(module, root, js)
        need(code == 1 and "absent.py 0 0 0.000000" in text, "missing python path was skipped")

        outside = root.parent / "outside-coverage.py"
        outside.write_text("x = 1\n")
        inner = root / "inner"
        inner.mkdir()
        payload = {"files": {str(outside): {"executed_lines": [1], "missing_lines": []}}}
        js = inner / "cov.json"
        js.write_text(json.dumps(payload))
        code, text = run(module, inner, js)
        need(code == 1 and "0 0 0.000000" in text, "outside path was skipped")
        present = inner / "seen.py"
        present.write_text("x = 1\n")
        payload = {"files": {str(present): {"executed_lines": [1], "missing_lines": []}}}
        js.write_text(json.dumps(payload))
        code, text = run(module, inner, js)
        need(code == 0 and "seen.py 1 1 1.000000" in text, "present python path was not scored")
        blank = root / "blank.cpp"
        blank.write_text("int z;\n")
        blank_gcov = root / "blank.cpp.gcov"
        blank_gcov.write_text("        -:    0:Source:blank.cpp\nnot a count\n        -:    1:int z;\n")
        code, text = run(module, root, blank_gcov)
        need("blank.cpp 0 0 1.000000" in text, "empty executable gcov was not accepted")
        empty_gcov = root / "empty.cpp.gcov"
        empty_gcov.write_text("        -:    0:Source:empty.cpp\n")
        (root / "empty.cpp").write_text("int only;\n")
        code, text = run(module, root, empty_gcov)
        need(code == 1, "gcov without source text was accepted")
        long_src = root / "long.cpp"
        long_src.write_text("int real;\n")
        long_gcov = root / "long.cpp.gcov"
        long_gcov.write_text("        -:    0:Source:long.cpp\n        1:    9:int real;\n")
        code, text = run(module, root, long_gcov)
        need(code == 1, "a line past the file was accepted")
        hidden = root / "bbdir"
        hidden.mkdir()
        (nested / "widget.cpp").write_text("int other;\n")
        (hidden / "widget.cpp").write_text("int real;\n")
        code, text = run(module, root, gcov)
        need(code == 1 and "0 0 0.000000" in text, "a build-directory copy was scored")
        comment_only = inner / "comment.py"
        comment_only.write_text("# only\n")
        payload = {"files": {str(comment_only): {"executed_lines": [], "missing_lines": []}}}
        js.write_text(json.dumps(payload))
        code, text = run(module, inner, js)
        need("comment.py 0 0 1.000000" in text, "comment-only python was not an empty set")
        note = root / "note.txt"
        note.write_text("x\n")
        code, text = run(module, root, note)
        need(code == 1 and "0 0 0.000000" in text, "plain text was accepted")
        outside.unlink()
        repo = root / "repo"
        repo.mkdir()
        env = os.environ.copy()
        env.update({
            "GIT_AUTHOR_NAME": "Rohit Goswami",
            "GIT_AUTHOR_EMAIL": "rgoswami@ieee.org",
            "GIT_COMMITTER_NAME": "Rohit Goswami",
            "GIT_COMMITTER_EMAIL": "rgoswami@ieee.org",
        })
        subprocess.run(["git", "init"], cwd=repo, check=True, capture_output=True)
        keep = repo / "keep.py"
        keep.write_text("x = 1\n")
        subprocess.run(["git", "add", "keep.py"], cwd=repo, check=True)
        first = subprocess.run(
            ["git", "-c", "core.hooksPath=/dev/null", "commit", "-m", "base"],
            cwd=repo, check=True, capture_output=True, text=True, env=env,
        )
        base = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
        (repo / "gap.py").write_text("y = 1\n")
        keep.write_text("x = 2\n")
        subprocess.run(["git", "add", "keep.py", "gap.py"], cwd=repo, check=True)
        subprocess.run(
            ["git", "-c", "core.hooksPath=/dev/null", "commit", "-m", "change"],
            cwd=repo, check=True, capture_output=True, env=env,
        )
        present = {"files": {str(keep): {"executed_lines": [1], "missing_lines": []}}}
        js = repo / "cov.json"
        js.write_text(json.dumps(present))
        saved = os.environ.get("DSEAMS_COVERAGE_BASE")
        os.environ["DSEAMS_COVERAGE_BASE"] = base
        try:
            code, text = run(module, repo, js)
            need(code == 1 and "gap.py 0 0 0.000000" in text, "a changed file with no row passed")
            present["files"][str(repo / "gap.py")] = {"executed_lines": [1], "missing_lines": []}
            js.write_text(json.dumps(present))
            code, text = run(module, repo, js)
            need(code == 0 and "gap.py 1 1 1.000000" in text, "a changed file with a row failed")
        finally:
            # saved was read after DSEAMS_COVERAGE_BASE was set to the
            # anchor commit, and nothing clears it before that read.
            os.environ["DSEAMS_COVERAGE_BASE"] = saved
        del first
        bare = Path(tempfile.mkdtemp())
        seen_py = bare / "seen.py"
        seen_py.write_text("x = 1\n")
        bare_js = bare / "cov.json"
        bare_js.write_text(json.dumps({"files": {str(seen_py): {"executed_lines": [1], "missing_lines": []}}}))
        os.environ.pop("DSEAMS_COVERAGE_BASE", None)
        code, text = run(module, bare, bare_js)
        need(code == 1 and "changed-files 7f310225 0 0 0.000000" in text, "a non-repository root passed")
        os.environ["DSEAMS_COVERAGE_BASE"] = "0123456789abcdef0123456789abcdef01234567"
        code, text = run(module, root, js)
        need(
            code == 1 and "changed-files 0123456789abcdef0123456789abcdef01234567 0 0 0.000000" in text,
            "an unknown base passed",
        )
        os.environ["DSEAMS_COVERAGE_BASE"] = ""
        code, text = run(module, root, js)
        need(code == 1 and "changed-files 0 0 0.000000" in text, "an empty base passed")
        real_git = subprocess.check_output(["sh", "-c", "command -v git"], text=True).strip()
        gitbin = root / "gitbin"
        gitbin.mkdir()
        wrapper = gitbin / "git"
        wrapper.write_text(
            "#!/bin/sh\n"
            "if [ \"$1\" = diff ]; then\n"
            "  exit 2\n"
            "fi\n"
            "exec %s \"$@\"\n" % real_git
        )
        wrapper.chmod(0o755)
        saved_path = os.environ["PATH"]
        os.environ["PATH"] = str(gitbin) + os.pathsep + saved_path
        os.environ["DSEAMS_COVERAGE_BASE"] = anchor
        try:
            code, text = run(module, root, js)
            need(code == 1 and ("changed-files " + anchor + " 0 0 0.000000") in text, "a git diff error passed")
        finally:
            os.environ["PATH"] = saved_path
        if saved_base is None:
            os.environ.pop("DSEAMS_COVERAGE_BASE", None)
        else:
            os.environ["DSEAMS_COVERAGE_BASE"] = saved_base
    print("coverage scorer ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
