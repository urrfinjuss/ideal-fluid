#include "ffluid.h"

void ffluid_disk() {
  printf("Module Loaded:\tDummy Disk\t\t");
  printf("Default precision:\tdouble\n");
}

void ffluid_setup_grid(data_ptr in) {
  unsigned long N = in->N;
  long_double_t l = in->l;
  long_double_t l2 = l*l;
  long_double_t q;
  for (unsigned long j = 0; j < N; j++) {
    q = PI*((2.0L*j/N) - 1.0L) - in->q0;
    in->q[j] = q + in->q0;
    in->u[j] = in->u0+2.L*atan2l(l*sinl(0.5L*q),cosl(0.5L*q));
    in->du[j] = -1.0IL*cexpl(-1.0IL*(q - in->q0));  // d\xi/ du
  }
}
