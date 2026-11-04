# Figure products the paper reads and git does not store.
# python3 paper/tex/figures/stats.py writes the plot tables and, when
# ring-head2head.json is present, the ring table and its macros.

PLOTS = [
    "plot-chill-dseams-chill.dat",
    "plot-chill-ovito-chill.dat",
    "plot-chill-palmer-chill.dat",
    "plot-churn.dat",
    "plot-nucleation.dat",
    "plot-ovrings-confined-dseams-rings.dat",
    "plot-ovrings-confined-ovito-rings.dat",
    "plot-ovrings-ic-dseams-rings.dat",
    "plot-ovrings-ic-ovito-rings.dat",
    "plot-ovrings-ih-dseams-rings.dat",
    "plot-ovrings-ih-ovito-rings.dat",
    "plot-ovrings-liquid-dseams-rings.dat",
    "plot-ovrings-liquid-ovito-rings.dat",
    "plot-q6-dseams-q6.dat",
    "plot-q6-freud-q6.dat",
    "plot-q6-pyscal-q6.dat",
    "plot-rings-dseams-rings.dat",
    "plot-rings-matscipy-rings-kernel.dat",
    "plot-seeded-dseams-seeded.dat",
    "plot-seeded-ovito-cna.dat",
    "plot-seeded-ovito-diamond.dat",
    "plot-seeded-ovito-ptm.dat",
    "plot-seeded-pyscal-diamond.dat",
    "plot-versions-one-ddc.dat",
    "plot-versions-one-hc.dat",
    "plot-versions-one-path.dat",
    "plot-versions-one-ring.dat",
    "plot-versions-two-ddc.dat",
    "plot-versions-two-hc.dat",
    "plot-versions-two-path.dat",
    "plot-versions-two-ring.dat",
]


rule all:
    input:
        expand("paper/tex/figures/data/{name}", name=PLOTS),
        "paper/tex/figures/data/plot-serial-rings.dat",
        "paper/tex/figures/data/ringperf_macros.tex",
        "paper/tex/figures/data/tab-ring-perf.tex",
        "paper/tex/figures/data/h2h/input-manifest.json",
        "paper/tex/figures/data/nucleation-host.txt",
        "paper/tex/figures/xyzrender_xyz/si_hydrate.xyz",
        "paper/tex/figures/xyzrender_xyz/sii_hydrate.xyz",


rule machine_records:
    output:
        "paper/tex/figures/data/host.txt",
        "paper/tex/figures/data/churn-host.txt",
        "paper/tex/figures/data/h2h/environment.json",
    shell:
        "python3 scripts/write_machine_records.py"


rule nucleation_host:
    output:
        "paper/tex/figures/data/nucleation-host.txt",
    shell:
        "python3 paper/tex/figures/generate_figure_data.py host"


rule input_manifest:
    output:
        "paper/tex/figures/data/h2h/input-manifest.json",
    shell:
        "python3 paper/tex/figures/generate_figure_data.py manifest"


rule hydrate_xyz:
    input:
        "input/traj/genice_sI.lammpstrj",
        "input/traj/genice_sII.lammpstrj",
    output:
        "paper/tex/figures/xyzrender_xyz/si_hydrate.xyz",
        "paper/tex/figures/xyzrender_xyz/sii_hydrate.xyz",
    shell:
        "python3 paper/tex/figures/generate_figure_data.py hydrate"


rule plots:
    input:
        "paper/tex/figures/stats.py",
        "paper/tex/figures/data/nucleation-host.txt",
        "paper/tex/figures/data/host.txt",
        "paper/tex/figures/data/churn-host.txt",
        "paper/tex/figures/data/h2h/environment.json",
    output:
        expand("paper/tex/figures/data/{name}", name=PLOTS),
    shell:
        "python3 paper/tex/figures/stats.py"


rule ring_panel:
    input:
        "paper/tex/figures/data/ring-head2head.json",
        "paper/tex/figures/stats.py",
    output:
        "paper/tex/figures/data/plot-serial-rings.dat",
        "paper/tex/figures/data/ringperf_macros.tex",
        "paper/tex/figures/data/tab-ring-perf.tex",
    shell:
        "python3 paper/tex/figures/stats.py"
