#!/usr/bin/env python3
"""Write supp_xref.tex and main_xref.tex from the .aux files.

A label's number is the first brace group of \\newlabel. rg_supplement.aux
feeds the s-prefix (\\sref in the main text). rg_main.aux feeds the
m-prefix (\\mref in the supplement).
"""
import re
import sys
from pathlib import Path

LABEL = re.compile(r"\\newlabel\{([^}]+)\}\{\{([^{}]*)\}")

def defs(aux: Path, prefix: str) -> str:
    lines = [
        "%% Written by scripts/mkxref.py from %s." % aux.name,
        "% Do not edit.",
        "% snapper:off",
        "\\makeatletter",
    ]
    text = aux.read_text(errors="replace")
    for name, number in sorted(dict(LABEL.findall(text)).items()):
        lines.append("\\xr@def{%s}{%s}{%s}" % (prefix, name, number))
    lines.append("\\makeatother")
    lines.append("% snapper:on")
    lines.append("")
    return "\n".join(lines)

def main() -> int:
    if len(sys.argv) != 4:
        print("usage: mkxref.py rg_supplement.aux rg_main.aux outdir", file=sys.stderr)
        return 2
    supp, main, out = map(Path, sys.argv[1:])
    (out / "supp_xref.tex").write_text(defs(supp, "s"))
    (out / "main_xref.tex").write_text(defs(main, "m"))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
