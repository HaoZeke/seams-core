.PHONY: check

# Fail when paper/tex/brine_science.tex and the eight gate CSVs disagree.
check:
	python3 input/brine/fit_velocity.py --check
	python3 -m pytest paper/data/brine-frames/test_brine.py -q
