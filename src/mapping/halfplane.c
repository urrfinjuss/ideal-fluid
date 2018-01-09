#include "ffluid.h"

void ffluid_halfplane() {
  printf("Module Loaded:\tDummy Half-Plane\t");
  printf("Default precision:\tdouble\n");
}

void ffluid_halfplane_setup_grid(data_ptr in) {
  unsigned long N = in->N;
  long_double_t l = in->l;
  long_double_t q;
  for (unsigned long j = 0; j < N; j++) {
    q = PI*((2.L*j/N) - 1.L) - in->q0;
    in->q[j] = q + in->q0; 
    in->u[j] = in->u0+2.L*atan2l(l*sinl(0.5L*q),cosl(0.5L*q));
    in->du[j] = 2.0L*l/(1.L + powl(l,2) + (1.L - powl(l,2))*cosl(q));
  }
}
