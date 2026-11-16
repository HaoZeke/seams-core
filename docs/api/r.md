# R

d-SEAMS 2.14.0. The R package is `bindings/r/dseams`. It calls
`libseams_lang` through cpp11. The exported frames are `read_frame`,
`nframes`, `positions`, `chill_plus`, `cages`, `primitive_rings`,
`signature_cages`, `symmetry`, `steinhardt`, `q4_q6`, `jump_tau90`, and
`rings_chordless`.
`print` on a frame dispatches to `print.dseams_frame`.
The methods are `positions.dseams_frame`, `chill_plus.dseams_frame`,
`cages.dseams_frame`, `primitive_rings.dseams_frame`,
`signature_cages.dseams_frame`, `symmetry.dseams_frame`,
`steinhardt.dseams_frame`, `q4_q6.dseams_frame`, and
`rings_chordless.dseams_frame`.

`docs/examples/core_counts.R` checks Ice Ih (128 hexagonal, 256
six-rings) and sII (136 clathrate, 144 five-rings, 16 six-rings).
`meson test docs_core_counts` runs it.

```{literalinclude} ../../../docs/examples/core_counts.R
:language: r
```

Build the library, then install the package against it:

```bash
meson setup build -Dwith_tests=true
meson compile -C build
DSEAMS_ROOT=$PWD R CMD INSTALL bindings/r/dseams
```

`configure` looks for `libseams_lang` under `build`, `bbdir`, `bbdir-py`, and
`build-eng`. For another build directory set `SEAMS_LIBDIR` to its
`bindings` directory; `libseams` and `seams_version.h` are then taken from the
`src` directory next to it.

`meson test r_testthat` installs the package from the build tree into
`r-test-lib` and runs the testthat file. CI runs the same file when `Rscript`
is on `PATH` and `DSEAMS_ROOT` points at this tree.
`docs/examples/minimal_r.R` reads the GenIce sI
cell and calls `steinhardt`, `chill_plus`, `primitive_rings`, and
`rings_chordless`.
`meson test docs_examples` runs it.

```{literalinclude} ../../../docs/examples/minimal_r.R
:language: r
```
`meson test docs_api_r` installs the package against `libseams_lang` and
runs the script below.

```{literalinclude} ../../../docs/examples/r_api.R
:language: r
```

`seams_lang_fixtures` is the C++ driver of `seams_lang.h`, not the R package.

```bash
DSEAMS_ROOT=$PWD meson test -C build language_fronts --print-errorlogs
```

```{literalinclude} ../../../bindings/r/dseams/tests/testthat/test-ice.R
:language: r
:lines: 9-36,77-91
```

`rings_chordless(frame, size = 6L, cutoff = 3.5, type = 0L, graph = "cutoff")`
lists the chordless cycles of one size, 3 through 12: cycles with no bond
between two of their atoms that are not next to each other on the cycle. It
returns an integer matrix with one row per cycle and `size` columns of 1-based
atom rows, so `positions(frame)[ring, ]` holds the coordinates of one cycle.
Each row starts at its smallest atom and the rows are sorted. The `per_atom`
attribute counts the cycles through each atom. `seams_rings_chordless(sys,
cutoff, size, type, graph)` takes a `seams_system_new` handle and returns
`list(rings, per_atom)` with the 0-based rows of the C call. `r_testthat`
compares both with the rows `chordless_ref` writes from the C API.

```{literalinclude} ../../../bindings/r/dseams/tests/testthat/test-ice.R
:language: r
:lines: 93-132
```

<!-- API_SYMBOLS_BEGIN -->

## Symbol reference

### cages

Purpose: Cage counts for a dump path.

Parameters: path; frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom; frame index is dimensionless

Example: `docs/examples/minimal_r.R`

### cages.dseams_frame

Purpose: Cage counts for a frame object.

Parameters: frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_r.R`

### chill_plus

Purpose: CHILL+ labels for a dump path.

Parameters: path; frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom

Example: `docs/examples/minimal_r.R`

### chill_plus.dseams_frame

Purpose: CHILL+ labels for a frame object.

Parameters: frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom; labels are dimensionless codes

Example: `docs/examples/minimal_r.R`

### jump_tau90

Purpose: Rotation time of H-H vectors between two frames.

Parameters: first; second; dt

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: dt and the result in the caller's time unit

Example: `docs/examples/minimal_r.R`

### nframes

Purpose: Number of frames in a trajectory path.

Parameters: path

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: dimensionless count

Example: `docs/examples/minimal_r.R`

### positions

Purpose: Coordinates of a dump frame.

Parameters: path; frame

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom

Example: `docs/examples/minimal_r.R`

### positions.dseams_frame

Purpose: Coordinates stored on a frame object.

Parameters: frame

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: angstrom

Example: `docs/examples/minimal_r.R`

### primitive_rings

Purpose: Ring histogram for a dump path.

Parameters: path; frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom

Example: `docs/examples/minimal_r.R`

### primitive_rings.dseams_frame

Purpose: Ring histogram for a frame object.

Parameters: frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_r.R`

### print.dseams_frame

Purpose: Print the atom count of a frame.

Parameters: frame

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: dimensionless count

Example: `docs/examples/minimal_r.R`

### q4_q6

Purpose: Local q4 and q6 for a dump path.

Parameters: path; frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_r.R`

### q4_q6.dseams_frame

Purpose: Local q4 and q6 for a frame object.

Parameters: frame; cutoff

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_r.R`

### read_frame

Purpose: Read one frame from a dump.

Parameters: path; frame; type

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: frame index is dimensionless; coordinates are angstrom

Example: `docs/examples/minimal_r.R`

### rings_chordless

Purpose: Chordless cycles of one size, one 1-based row per cycle.

Parameters: frame; size; cutoff; type; graph

Units: cutoff in angstrom; indices and counts are dimensionless

Example: `docs/examples/minimal_r.R`

### rings_chordless.dseams_frame

Purpose: Chordless cycles for a frame object.

Parameters: frame; size; cutoff; type; graph

Units: cutoff in angstrom; per_atom counts are dimensionless

Example: `docs/examples/minimal_r.R`

### signature_cages

Purpose: Closed cages for a dump path.

Parameters: path; frame; spec

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: spec is text; cutoff used by the call is in angstrom

Example: `docs/examples/minimal_r.R`

### signature_cages.dseams_frame

Purpose: Closed cages for a frame object.

Parameters: frame; spec

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: spec is text

Example: `docs/examples/minimal_r.R`

### steinhardt

Purpose: Mean local q_l for a dump path.

Parameters: path; frame; order_l

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: order parameter is dimensionless

Example: `docs/examples/minimal_r.R`

### steinhardt.dseams_frame

Purpose: Mean local q_l for a frame object.

Parameters: frame; order_l

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: order parameter is dimensionless

Example: `docs/examples/minimal_r.R`

### symmetry

Purpose: Shell labels for a dump path.

Parameters: path; frame

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: dimensionless label codes

Example: `docs/examples/minimal_r.R`

### symmetry.dseams_frame

Purpose: Shell labels for a frame object.

Parameters: frame

Returns: The data frame, vector, or count named in the purpose.

Ownership: R owns the returned vectors. A dseams_frame finalizer frees the C frame.

Errors: The binding stops with the text of seams_last_error when the status is not SEAMS_OK.

Units: dimensionless label codes

Example: `docs/examples/minimal_r.R`

### seams_system_new

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_system_read

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_system_count

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_chill

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_steinhardt

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_cages

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_signature_cages

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_neighbours_knn

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_rdf

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_running_cn

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_coordination

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_hbonds

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_density_axis

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_fingerprint

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_contact_pairs

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_domains

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_prisms

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_monolayer

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_rings_vertices

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`

### seams_site_table

Purpose: Call into libseams through the C shim.

Parameters: see include/seams.h

Units: angstrom where a length is involved; otherwise dimensionless

Example: `docs/examples/r_api.R`
<!-- API_SYMBOLS_END -->


## Structure still

```{figure} ../tutorials/molecules/r/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug r` writes the 2400 by 1800 PNG.
```

