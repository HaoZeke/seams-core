# Head-to-head watch. Frames are generated at run time and are not stored.
# d-SEAMS is bbdir-release/tests/bench_watch from meson -Dbuildtype=release.

rule related_watch_v3:
    input:
        "scripts/h2h_v3.py",
        "scripts/h2h_v2.py",
        "scripts/h2h_mdapy.py",
        "scripts/h2h_protocol.py",
        "bbdir-release/tests/bench_watch",
    output:
        "benchmarks/h2h-watch-v3.csv",
        "benchmarks/h2h-watch-v3.log",
    shell:
        "PY=\"${{H2H_PYTHON:-python3}}\"; "
        "\"$PY\" scripts/h2h_v3.py "
        "--bench bbdir-release/tests/bench_watch "
        "--python \"$PY\" "
        "--out benchmarks/h2h-watch-v3.csv "
        "--log benchmarks/h2h-watch-v3.log"
