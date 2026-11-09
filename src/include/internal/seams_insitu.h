#ifndef SEAMS_INSITU_H_
#define SEAMS_INSITU_H_

#ifdef __cplusplus
extern "C" {
#endif

/* In-situ frame analysis for a running molecular dynamics calculation.
 *
 * Coordinates are Angstrom. The box is the LAMMPS edge frame: xlo, ylo, zlo
 * is the edge origin, lx, ly, lz are the edge lengths, and xy, xz, yz are
 * the tilt factors. A zero tilt is orthorhombic.
 *
 * xyz, types and ids are borrowed for the call and are not retained.
 * The engine copies them into its point cloud and analyses that copy.
 * The neighbour list, the primitive rings, the published CHILL+ labels
 * and the cage flags are computed from those values. The hot path does
 * not write a file and does not read one back. debug_dump, off by
 * default, appends a LAMMPS dump after the analysis for inspection.
 *
 * gpu_resident means xyz is a device pointer. This build does not
 * dereference it and does not copy device memory. The host copy is left
 * for a later run.
 *
 * local_only means the arrays hold one rank's atoms. The adapter gathers
 * the full system with MPI_Gatherv onto the analysis rank and calls
 * again with gathered set. The core does not link that library. A call
 * that is still local returns SEAMS_INSITU_NEED_GATHER.
 *
 * CHILL+ is the published classifier. This header does not change it.
 */

#define SEAMS_INSITU_OK 0
#define SEAMS_INSITU_ERR 1
#define SEAMS_INSITU_SKIP 2
#define SEAMS_INSITU_NEED_GATHER 3
#define SEAMS_INSITU_GPU_DEFERRED 4
#define SEAMS_INSITU_ERR_BOX 5
#define SEAMS_INSITU_ERR_IO 6

#define SEAMS_INSITU_RING_MAX 12
#define SEAMS_INSITU_CHILL_N 9

typedef struct seams_insitu_frame {
  int step;
  int n;
  const double *xyz;
  const int *types;
  const int *ids;
  double xlo, ylo, zlo;
  double lx, ly, lz;
  double xy, xz, yz;
  int local_only;
  int gathered;
  int gpu_resident;
  int gpu_device;
  int rank;
  int nranks;
} seams_insitu_frame;

typedef struct seams_insitu_cadence {
  int every;
  int first;
  int last;
} seams_insitu_cadence;

typedef struct seams_insitu_options {
  double cutoff;
  int type;
  int ring_depth;
  int write_labels;
  seams_insitu_cadence cadence;
  /* 0 analyses memory only. 1 also appends a LAMMPS dump when a path
   * is passed to seams_insitu_observe. */
  int debug_dump;
} seams_insitu_options;

typedef struct seams_insitu_result {
  int status;
  int step;
  int n;
  int n_rings;
  int ring_count[SEAMS_INSITU_RING_MAX + 1];
  int chill_count[SEAMS_INSITU_CHILL_N];
  int n_hc;
  int n_ddc;
  int n_both;
  int n_water;
  int n_water_atoms;
  const int *chill_labels;
  const int *cage_codes;
  const char *text;
  int text_len;
} seams_insitu_result;

typedef int (*seams_insitu_sink)(const seams_insitu_result *result, void *user);

typedef struct seams_insitu_engine seams_insitu_engine;

void seams_insitu_options_default(seams_insitu_options *opt);

/* Null when the options are rejected. Null options take the defaults:
 * cutoff 3.5 Angstrom, water type 1, ring depth 8, labels on,
 * every step from 0 with no last step, debug dump off. */
seams_insitu_engine *seams_insitu_create(const seams_insitu_options *opt);
void seams_insitu_destroy(seams_insitu_engine *eng);

int seams_insitu_set_sink(seams_insitu_engine *eng, seams_insitu_sink sink,
                          void *user);
int seams_insitu_open_record(seams_insitu_engine *eng, const char *path);
int seams_insitu_close_record(seams_insitu_engine *eng);

/* 1 when step falls on the cadence, 0 when it does not, -1 when the
 * cadence itself is rejected. */
int seams_insitu_due(const seams_insitu_cadence *cadence, int step);

int seams_insitu_step(seams_insitu_engine *eng, const seams_insitu_frame *frame);

/* The next analysed frame uses this neighbour list instead of building
 * one. n rows, counts[i] long, concatenated in flat. Each row starts
 * with that atom's id. A wider host list is filtered to the cutoff.
 * n == 0 clears the override. */
int seams_insitu_set_neighbours(seams_insitu_engine *eng, int n,
                                const int *counts, const int *flat);

/* Analyses the borrowed frame. When debug_dump is set and dump_path is
 * not empty, the same frame is appended to that file after the analysis.
 * The dump is not read back. A null or empty path never writes. */
int seams_insitu_observe(seams_insitu_engine *eng, const seams_insitu_frame *frame,
                         const char *dump_path);

int seams_insitu_append_dump(const seams_insitu_frame *frame, const char *path);

/* Analyses every frame of a dump written by seams_insitu_append_dump or
 * seams_insitu_observe. The record file receives the same text the sink
 * would have written. The sink is not called. */
int seams_insitu_replay(seams_insitu_engine *eng, const char *dump_path,
                        const char *record_path);

int seams_insitu_status(const seams_insitu_engine *eng);
const seams_insitu_result *seams_insitu_last(const seams_insitu_engine *eng);

int seams_insitu_main(int argc, char **argv);

/* Interactive molecular dynamics client. box9 is the angstrom edge box.
 * max_frames <= 0 reads until the server closes. Positions arrive in angstrom. */
int seams_insitu_imd_session(const char *host, int port, const char *dump,
                             const char *record, const double *box9, int max_frames);
int seams_insitu_imd_main(int argc, char **argv);

int seams_insitu_frame_bytes(void);
int seams_insitu_options_bytes(void);
int seams_insitu_result_bytes(void);

/* Keep atoms whose mask bit is set. Null types or ids become type 1 and
 * ids 1..nlocal. Returns the kept count, or -1 when the arguments are
 * rejected. */
int seams_insitu_pack_owned(int nlocal, const int *mask, int groupbit,
                            const double *xyz, const int *types, const int *ids,
                            double *xyz_out, int *types_out, int *ids_out);

/* Exclusive prefix of counts. total receives the sum. */
int seams_insitu_gather_displs(int nranks, const int *counts, int *displs,
                               int *total);

/* Nanometre coordinates to Angstrom. */
void seams_insitu_scale_nm(const double *xyz_nm, int n, double *xyz_a);

/* Row-major lattice vectors in the reduced frame: a along x, b in the
 * xy plane. Origin is zero. Returns 0 on success. */
int seams_insitu_box_from_rows(const double rows[9], double *xlo, double *ylo,
                               double *zlo, double *lx, double *ly, double *lz,
                               double *xy, double *xz, double *yz);

/* hc and ddc are 0 or 1 per atom. codes are 0 water, 1 hexagonal cage,
 * 2 double-diamond cage, 3 both. Atoms whose type is not water_type
 * stay 0 and are left out of the four counters. */
void seams_insitu_tally_cages(int n, const int *types, int water_type,
                              const int *hc, const int *ddc, int *codes,
                              int *n_hc, int *n_ddc, int *n_both, int *n_water);

#ifdef __cplusplus
}
#endif

#endif
