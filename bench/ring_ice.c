/*
** SQQ's ice-I rule (sqq/core/sqq_py/ice.py, SQQ v0.5.8) on d-SEAMS rings,
** for the SQQ head-to-head in scripts/h2h_sqq_versus.py. Not part of the
** public API; it uses only seams.h.
**
** A molecule with four neighbours on at least one six-ring is interfacial.
** It is ice when it sits on at least min_six six-rings and every neighbour
** also has four neighbours. Every other molecule is liquid.
**
** definition 0 counts the primitive six-rings of seams_rings_find.
** definition 1 counts the chordless six-cycles of the same cutoff graph,
** which is the ring set SQQ's ice stage counts.
*/

#include "seams.h"

#include <stdlib.h>
#include <string.h>

enum { H2H_PRIMITIVE = 0, H2H_CHORDLESS = 1 };
enum { H2H_LIQUID = 0, H2H_INTERFACIAL = 1, H2H_ICE = 2 };

typedef struct {
  const int *off;
  const int *cols;
  int *six;
  int *rings;
  size_t cap;
  size_t *next;
} walk_t;

static void record(const walk_t *w, const int *path) {
  for (int k = 0; k < 6; ++k) {
#pragma omp atomic
    w->six[path[k]] += 1;
  }
  if (w->next != NULL) {
    size_t slot;
#pragma omp atomic capture
    slot = (*w->next)++;
    if (w->rings != NULL && slot < w->cap) {
      memcpy(w->rings + 6 * slot, path, 6 * sizeof(int));
    }
  }
}

/* mark[x] is one more than the position of x on the path, or 0. Each
** cycle is found once: from its lowest vertex, with path[1] < path[5]. */
static void extend(const walk_t *w, unsigned char *mark, int *path, int k) {
  const int s = path[0];
  const int last = path[k - 1];
  for (int a = w->off[last]; a < w->off[last + 1]; ++a) {
    const int v = w->cols[a];
    if (v <= s || mark[v] != 0) {
      continue;
    }
    int closes = 0;
    int chord = 0;
    for (int b = w->off[v]; b < w->off[v + 1] && !chord; ++b) {
      const int pos = mark[w->cols[b]];
      if (pos == 0 || pos == k) {
        continue;
      }
      if (pos == 1) {
        closes = 1;
      } else {
        chord = 1;
      }
    }
    if (chord) {
      continue;
    }
    if (k == 5) {
      if (closes && path[1] < v) {
        path[5] = v;
        record(w, path);
      }
      continue;
    }
    if (closes) {
      continue;
    }
    mark[v] = (unsigned char)(k + 1);
    path[k] = v;
    extend(w, mark, path, k + 1);
    mark[v] = 0;
  }
}

static int primitive_six(seams_system *sys, double cutoff, const walk_t *w) {
  seams_rings *found = NULL;
  int status = seams_rings_find(sys, cutoff, 0, SEAMS_GRAPH_CUTOFF, 6, &found);
  if (status != SEAMS_OK) {
    return status;
  }
  const int count = seams_rings_count(found);
  int *sizes = malloc((size_t)(count > 0 ? count : 1) * sizeof(int));
  if (sizes == NULL) {
    seams_rings_destroy(found);
    return SEAMS_ERR_NOMEM;
  }
  if (count > 0) {
    status = seams_rings_labels(found, sizes, NULL);
  }
  for (int k = 0; k < count && status == SEAMS_OK; ++k) {
    if (sizes[k] != 6) {
      continue;
    }
    int ring[6];
    size_t got = 0;
    status = seams_rings_vertices(found, k, ring, 6, &got);
    if (status == SEAMS_OK) {
      record(w, ring);
    }
  }
  free(sizes);
  seams_rings_destroy(found);
  return status;
}

static int chordless_six(int n, const walk_t *w) {
  int failed = 0;
#pragma omp parallel num_threads(seams_get_num_threads())
  {
    unsigned char *mark = calloc((size_t)(n > 0 ? n : 1), 1);
    if (mark == NULL) {
#pragma omp atomic write
      failed = 1;
    }
#pragma omp for schedule(dynamic, 256)
    for (int s = 0; s < n; ++s) {
      if (mark == NULL) {
        continue;
      }
      int path[6];
      path[0] = s;
      mark[s] = 1;
      extend(w, mark, path, 1);
      mark[s] = 0;
    }
    free(mark);
  }
  return failed ? SEAMS_ERR_NOMEM : SEAMS_OK;
}

/*
** labels and six take n ints: 0 liquid, 1 interfacial, 2 ice, and the
** number of counted six-rings through each molecule. rings, when not
** NULL, receives up to rings_cap rings of six vertex indices. n_rings,
** when not NULL, receives the number of counted six-rings; a short rings
** buffer returns SEAMS_ERR_RANGE after the labels are written.
*/
int seams_h2h_ring_ice(seams_system *sys, double cutoff, int definition, int min_six,
                       int *labels, int *six, size_t n, int *rings, size_t rings_cap,
                       size_t *n_rings) {
  if (sys == NULL || labels == NULL || six == NULL) {
    return SEAMS_ERR_NULL;
  }
  if (definition != H2H_PRIMITIVE && definition != H2H_CHORDLESS) {
    return SEAMS_ERR_ARG;
  }
  seams_neighbours *nb = NULL;
  int status = seams_neighbours_build(sys, cutoff, 0, SEAMS_GRAPH_CUTOFF, &nb);
  if (status != SEAMS_OK) {
    return status;
  }
  size_t n_off = 0;
  size_t nnz = 0;
  int *off = NULL;
  int *cols = NULL;
  status = seams_neighbours_copy(nb, NULL, 0, &n_off, NULL, 0, &nnz);
  if (status == SEAMS_OK && n_off != n + 1) {
    status = SEAMS_ERR_ARG;
  }
  if (status == SEAMS_OK) {
    off = malloc(n_off * sizeof(int));
    cols = malloc((nnz > 0 ? nnz : 1) * sizeof(int));
    status = off == NULL || cols == NULL
                 ? SEAMS_ERR_NOMEM
                 : seams_neighbours_copy(nb, off, n_off, &n_off, cols, nnz, &nnz);
  }
  seams_neighbours_destroy(nb);
  size_t next = 0;
  if (status == SEAMS_OK) {
    memset(six, 0, n * sizeof(int));
    const walk_t w = {off, cols, six, rings, rings_cap,
                      rings != NULL || n_rings != NULL ? &next : NULL};
    status = definition == H2H_PRIMITIVE ? primitive_six(sys, cutoff, &w)
                                         : chordless_six((int)n, &w);
  }
  if (status == SEAMS_OK) {
    const int count = (int)n;
#pragma omp parallel for num_threads(seams_get_num_threads()) schedule(static)
    for (int i = 0; i < count; ++i) {
      int label = H2H_LIQUID;
      if (off[i + 1] - off[i] == 4 && six[i] > 0) {
        label = H2H_INTERFACIAL;
        if (six[i] >= min_six) {
          int ordered = 1;
          for (int a = off[i]; a < off[i + 1]; ++a) {
            const int j = cols[a];
            ordered &= off[j + 1] - off[j] == 4;
          }
          if (ordered) {
            label = H2H_ICE;
          }
        }
      }
      labels[i] = label;
    }
    if (n_rings != NULL) {
      *n_rings = next;
    }
    if (rings != NULL && next > rings_cap) {
      status = SEAMS_ERR_RANGE;
    }
  }
  free(off);
  free(cols);
  return status;
}
