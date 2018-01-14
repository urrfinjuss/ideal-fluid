#include "ffluid.h"

void ffluid_disk() {
  printf("Module Loaded:\tDummy Disk\t\t");
  printf("Default precision:\tdouble\n");
}

void ffluid_setup_grid(data_ptr in) {
  unsigned long N = in->N;
  long_double_t l = in->l;
  long_double_t l2 = l*l;
  long_double_t q0 = in->q0;
  long_double_t q;
  for (unsigned long j = 0; j < N; j++) {
    q = PI*((2.0L*j/N) - 1.0L) - q0;
    in->q[j] = q + q0; 
    in->u[j] = 1.0L + 1.0IL*l*tanl(0.5L*q);
    in->u[j] = (1.0L - 1.0IL*l*tanl(0.5L*q))*cexpl(-1.0IL*in->u0)/in->u[j];
    in->du[j] = 1.0L - l2 + (1.0L + l2)*cosl(q) + 2.0IL*l*sinl(q);
    in->du[j] = -2.0IL*l*cexpl(-1.0IL*in->u0)/in->du[j];
  }
}
