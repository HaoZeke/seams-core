#ifndef SEAMS_C_API_H_
#define SEAMS_C_API_H_

#ifdef __cplusplus
extern "C" {
#endif

/* C ABI for native callers and for an emscripten/WASM build.
 * xyz is n*3, box is ox,oy,oz,lx,ly,lz (ortho). labels[i] is
 * 0 cubic, 1 hexagonal, 2 water, 3 interfacial, 4 clathrate,
 * 5 interClathrate, 6 unclassified, 7 reCubic, 8 reHex.
 * Returns 0 on success.
 * seams_chill_plus bonds every pair closer than 3.5 A (the cutoff list).
 * seams_chill_plus_graph takes the cutoff and graph: "cutoff" (also NULL
 * or "") or "mutual", the mutual four nearest inside cutoff + 1.5. */
int seams_chill_plus(const double *xyz, int n, const double *box, int *labels);
int seams_chill_plus_graph(const double *xyz, int n, const double *box,
                           double cutoff, const char *graph, int *labels);
/* rule is "chill+" (NULL or "") or "ovito-chill+". A YAML or Lua
 * document with a rule field is accepted. Returns 1 on a bad rule. */
int seams_chill_plus_rule(const double *xyz, int n, const double *box,
                          double cutoff, const char *graph, const char *rule,
                          int *labels);

/* box9 is xlo, ylo, zlo, lx, ly, lz, xy, xz, yz. lx, ly, lz are the
 * edge lengths xhi-xlo, not the bound spans, and xlo, ylo, zlo are the
 * edge origin. A zero tilt matches seams_chill_plus_graph. */
int seams_chill_plus_graph_lammps(const double *xyz, int n,
                                  const double *box9, double cutoff,
                                  const char *graph, int *labels);

/* Mixed frames. types and mols may be NULL (every atom is then one
 * species, mol id = index + 1). oxygen_type > 0 keeps that LAMMPS type.
 * mol_mode 0 ignores molecules, 1 keeps the water oxygen of each
 * multi-atom molecule, 2 keeps mol_ids[0..n_mol_ids). A mixed frame
 * with mol_mode 0 and oxygen_type 0 is an error. Labels of atoms that
 * are not water oxygens are unclassified (6). ion_ice, ion_liquid and
 * ion_front may be NULL. ion_front is the nearest ice oxygen, or -1.
 * seams_c_error() is the message after a non-zero return.
 * On a single-species cloud the labels match seams_chill_plus_graph_lammps. */
int seams_chill_plus_species(const double *xyz, int n, const double *box9,
                             const int *types, const int *mols,
                             int oxygen_type, int mol_mode, const int *mol_ids,
                             int n_mol_ids, double cutoff, const char *graph,
                             int *labels, int *ion_ice, int *ion_liquid,
                             double *ion_front);

/* degree[i] is the oxygen-oxygen neighbour count inside cutoff, or -1
 * when atom i is not a selected water oxygen. */
int seams_oxygen_cutoff_degree(const double *xyz, int n, const double *box9,
                               const int *types, const int *mols,
                               int oxygen_type, int mol_mode,
                               const int *mol_ids, int n_mol_ids, double cutoff,
                               int *degree);

const char *seams_c_error(void);

#ifdef __cplusplus
}
#endif

#endif
