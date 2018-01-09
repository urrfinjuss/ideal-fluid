#include "ffluid.h"

void ffluid_disk() {
  printf("Module Loaded:\tDummy Disk\t\t");
  printf("Default precision:\tdouble\n");
}

void ffluid_disc_setup_grid(data_ptr in) {
  unsigned long N = in->N;
  long_double_t l = in->l;
  long_double_t q;
  for (unsigned long j = 0; j < N; j++) {
    q = PI*((2.L*j/N) - 1.L);
    in->q[j] = q; 
    in->u[j] = cexpl(1.I*q);
    in->du[j] = 1.L;
  }
}
