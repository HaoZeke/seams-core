# Head-to-head watch. Frames are generated at run time and are not stored.
# d-SEAMS is the release build in bbdir-release (meson -Dbuildtype=release).

rule related_watch:
    input:
        "scripts/h2h_protocol.py",
        "scripts/h2h_watch.py",
        "scripts/h2h_mdapy.py",
        "scripts/h2h_fastatomstruct.py",
        "benchmarks/h2h_watch.py",
        "bbdir-release/src/libyodaLib.so",
    output:
        "benchmarks/h2h-watch.csv",
    shell:
        "PY=\"${{H2H_PYTHON:-python3}}\"; "
        "\"$PY\" scripts/h2h_watch.py "
        "--lib bbdir-release/src/libyodaLib.so "
        "--python \"$PY\" "
        "--out benchmarks/h2h-watch.csv"


rule related_watch_v2:
    input:
        "scripts/h2h_protocol.py",
        "scripts/h2h_watch.py",
        "scripts/h2h_v2.py",
        "bbdir-release/src/libyodaLib.so",
    output:
        "benchmarks/h2h-watch-v2.csv",
    shell:
        "PY=\"${{H2H_PYTHON:-python3}}\"; "
        "\"$PY\" scripts/h2h_v2.py "
        "--lib bbdir-release/src/libyodaLib.so "
        "--python \"$PY\" "
        "--out benchmarks/h2h-watch-v2.csv"
