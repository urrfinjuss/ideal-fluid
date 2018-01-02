#include "ffluid.h"

void ffluid_halfplane() {
  printf("Module Loaded:\tDummy Half-Plane\t");
  printf("Default precision:\tdouble\n");
}

void ffluid_setup_grid(data_ptr in, grid_ptr out) {
  unsigned long N = in->N;
  long_double_t	q;
  if (out->Initialized) ffluid_reinit_grid(in, out);
  else ffluid_init_grid(in, out);
  for (unsigned long j = 0; j < N; j++) {
    q = PI*((2*j/N) - 1);
    out->q[j] = q; 
    out->u[j] = in->u0 + 2.*atan2(in->l*sin(0.5*(q - in->q0)), cos(0.5*(q - in->q0)));
  }
}
