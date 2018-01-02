#include "ffluid.h"


void ffluid_memory_module() {
  printf("Module ffluid/memory.h:\n");
  printf("memory.c:\tffluid_memory_module ffluid_init_data ffluid_init_grid\n");
}

void ffluid_init_data(data_ptr in, unsigned long N) {
  in->Q = malloc(N*sizeof(long_complex_t));
  in->V = malloc(N*sizeof(long_complex_t));
}

void ffluid_init_grid(data_ptr in, grid_ptr out) {
  out->N = in->N;
  out->q = malloc(in->N*sizeof(long_double_t));
  out->u = malloc(in->N*sizeof(long_double_t));
}

void ffluid_reinit_grid(data_ptr in, grid_ptr out) {
  free(out->q);
  free(out->u);
  out->N = in->N;
  out->q = malloc(in->N*sizeof(long_double_t));
  out->u = malloc(in->N*sizeof(long_double_t));
}
