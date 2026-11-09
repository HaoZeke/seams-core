#include "seams_insitu.h"

#include <stdio.h>
#include <string.h>

static int rejected_cutoff(void) {
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  opt.cutoff = 0.0;
  if (seams_insitu_create(&opt) != NULL) return 0;
  fprintf(stderr, "create failed\n");
  return 1;
}

static int rejected_record(void) {
  seams_insitu_options opt;
  seams_insitu_engine *eng;
  seams_insitu_options_default(&opt);
  eng = seams_insitu_create(&opt);
  if (eng == NULL) return 0;
  if (seams_insitu_open_record(eng, "/no/such/dir/census.txt") == 0) { seams_insitu_destroy(eng); return 0; }
  fprintf(stderr, "record failed\n");
  seams_insitu_destroy(eng);
  return 2;
}

static int rejected_observe(void) {
  seams_insitu_options opt;
  seams_insitu_frame frame;
  seams_insitu_engine *eng;
  double xyz[3] = {0.0, 0.0, 0.0};
  int types[1] = {1};
  int rc;
  seams_insitu_options_default(&opt);
  opt.debug_dump = 1;
  eng = seams_insitu_create(&opt);
  if (eng == NULL) return 0;
  memset(&frame, 0, sizeof frame);
  frame.n = 1;
  frame.xyz = xyz;
  frame.types = types;
  frame.gathered = 1;
  rc = seams_insitu_observe(eng, &frame, "frame.lammpstrj");
  seams_insitu_destroy(eng);
  if (rc == SEAMS_INSITU_OK) return 0;
  fprintf(stderr, "observe status %d\n", rc);
  return 3;
}

int main(void) {
  double xyz[12] = {0.0, 0.0, 0.0, 2.8, 0.0, 0.0, 0.0, 2.8, 0.0, 1.4, 1.4, 2.0};
  int types[4] = {1, 1, 1, 1};
  seams_insitu_options opt;
  seams_insitu_frame frame;
  seams_insitu_engine *eng;
  int rc;

  if (rejected_cutoff() != 1) return 1;
  if (rejected_record() != 2) return 2;
  if (rejected_observe() != 3) return 3;
  seams_insitu_options_default(&opt);
  opt.cutoff = 3.5;
  opt.type = 1;
  opt.debug_dump = 1;
  eng = seams_insitu_create(&opt);
  if (eng == NULL) return 1;
  if (seams_insitu_open_record(eng, "census.txt") != 0) return 2;
  memset(&frame, 0, sizeof frame);
  frame.step = 0;
  frame.n = 4;
  frame.xyz = xyz;
  frame.types = types;
  frame.lx = frame.ly = frame.lz = 20.0;
  frame.gathered = 1;
  frame.coord_scale = 1.0;
  rc = seams_insitu_observe(eng, &frame, "frame.lammpstrj");
  seams_insitu_close_record(eng);
  seams_insitu_destroy(eng);
  if (rc != SEAMS_INSITU_OK) return 3;
  printf("observe %d\n", rc);
  return 0;
}
