#include "seams_insitu.h"

#include <stdio.h>
#include <string.h>

int main(void) {
  double xyz[12] = {0.0, 0.0, 0.0, 2.8, 0.0, 0.0, 0.0, 2.8, 0.0, 1.4, 1.4, 2.0};
  int types[4] = {1, 1, 1, 1};
  seams_insitu_options opt;
  seams_insitu_frame frame;
  seams_insitu_engine *eng;
  int rc;

  seams_insitu_options_default(&opt);
  opt.cutoff = 3.5;
  opt.type = 1;
  opt.debug_dump = 1;
  eng = seams_insitu_create(&opt);
  if (eng == NULL) {
    fprintf(stderr, "create failed\n");
    return 1;
  }
  if (seams_insitu_open_record(eng, "census.txt") != 0) {
    fprintf(stderr, "record failed\n");
    return 2;
  }
  memset(&frame, 0, sizeof(frame));
  frame.step = 0;
  frame.n = 4;
  frame.xyz = xyz;
  frame.types = types;
  frame.lx = 20.0;
  frame.ly = 20.0;
  frame.lz = 20.0;
  frame.gathered = 1;
  frame.coord_scale = 1.0;
  rc = seams_insitu_observe(eng, &frame, "frame.lammpstrj");
  seams_insitu_close_record(eng);
  seams_insitu_destroy(eng);
  if (rc != SEAMS_INSITU_OK) {
    fprintf(stderr, "observe status %d\n", rc);
    return 3;
  }
  printf("observe %d\n", rc);
  return 0;
}
