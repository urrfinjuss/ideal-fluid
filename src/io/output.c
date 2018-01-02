#include "ffluid.h"

void ffluid_write_surface(data_ptr in, char *fname) {
  unsigned long N = in->N;
  __float128	time = in->time;
  char path[80], input_file[80];
  printf("%s\n", Control.data_path);

  for (unsigned long j = 0; j < N; j++) {
  }
}

void ffluid_append_to_log() {
}
