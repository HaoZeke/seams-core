/* Copyright (c) d-SEAMS developers */
/* seams_rings_chordless through the public C API, written as whitespace
** separated integers for the binding tests to compare against:
**   atoms n_rings size
**   n_rings rows of size 0-based atom indices
**   one per-atom cycle count per atom
** usage: chordless_ref path frame type_filter graph cutoff type size [out] */
#include "seams.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int fail(const char *what) {
  fprintf(stderr, "chordless_ref: %s: %s\n", what, seams_last_error());
  return 1;
}

int main(int argc, char **argv) {
  if (argc != 8 && argc != 9) {
    fprintf(stderr, "usage: chordless_ref path frame type_filter cutoff|mutual cutoff type size [out]\n");
    return 2;
  }
  const char *path = argv[1];
  const int frame = atoi(argv[2]);
  const int type_filter = atoi(argv[3]);
  const int graph = strcmp(argv[4], "mutual") == 0 ? SEAMS_GRAPH_MUTUAL : SEAMS_GRAPH_CUTOFF;
  const double cutoff = atof(argv[5]);
  const int type = atoi(argv[6]);
  const int size = atoi(argv[7]);
  FILE *out = stdout;
  if (argc == 9) {
    out = fopen(argv[8], "w");
    if (out == NULL) {
      perror(argv[8]);
      return 1;
    }
  }
  seams_system *sys = seams_system_create();
  if (sys == NULL) {
    return fail("system");
  }
  if (seams_system_read(sys, path, frame, type_filter) != SEAMS_OK) {
    return fail("read");
  }
  seams_neighbours *nb = NULL;
  if (seams_neighbours_build(sys, cutoff, type, graph, &nb) != SEAMS_OK) {
    return fail("graph");
  }
  const int n = seams_system_count(sys);
  int *per_atom = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
  size_t count = 0;
  if (seams_rings_chordless(nb, size, NULL, 0, &count, per_atom, (size_t)n) != SEAMS_OK) {
    return fail("count");
  }
  int *rings = malloc((count > 0 ? count : 1) * (size_t)size * sizeof(int));
  if (seams_rings_chordless(nb, size, rings, count * (size_t)size, &count, NULL, 0) != SEAMS_OK) {
    return fail("rings");
  }
  fprintf(out, "%d %zu %d\n", n, count, size);
  for (size_t r = 0; r < count; ++r) {
    for (int k = 0; k < size; ++k) {
      fprintf(out, k == 0 ? "%d" : " %d", rings[r * (size_t)size + (size_t)k]);
    }
    fputc('\n', out);
  }
  for (int i = 0; i < n; ++i) {
    fprintf(out, "%d\n", per_atom[i]);
  }
  free(rings);
  free(per_atom);
  seams_neighbours_destroy(nb);
  seams_system_destroy(sys);
  if (out != stdout && fclose(out) != 0) {
    perror(argv[8]);
    return 1;
  }
  return 0;
}
