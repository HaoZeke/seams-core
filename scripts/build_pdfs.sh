#!/bin/sh
# Build both PDFs. Exit 1 when supp_xref.tex or main_xref.tex is stale, or
# when a log reports an undefined reference, an undefined citation or an
# overfull box.
set -eu
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
tex=$root/paper/tex
tectonic=${TECTONIC:-tectonic}
cd "$tex"
# Tectonic 0.15 sees the elsarticle .bbl as changed on every pass and stops
# at its rerun limit, so the passes after BibTeX are fixed.
# RERUNS passes after the first settle the citations, the cross-references
# and the page numbers of both documents.
RERUNS=3
for doc in rg_supplement rg_main; do
  "$tectonic" -X compile --keep-intermediates --keep-logs -r "$RERUNS" \
    "$doc.tex"
done
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
python3 "$root/scripts/mkxref.py" rg_supplement.aux rg_main.aux "$tmp"
stale=0
for name in supp_xref.tex main_xref.tex; do
  if ! cmp -s "$name" "$tmp/$name"; then
    echo "stale: paper/tex/$name" >&2
    diff -u "$name" "$tmp/$name" >&2 || true
    stale=1
  fi
done
for doc in rg_supplement rg_main; do
  if grep -E "undefined|Overfull|Rerun to get" "$doc.log" >&2; then
    echo "unclean: paper/tex/$doc.log" >&2
    stale=1
  fi
done
exit "$stale"
