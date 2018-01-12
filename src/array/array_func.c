#include "ffluid.h"

void ffluid_arrayfunc_module() {
  printf("Module ffluid/array.h:\n");
  printf("array_func.c:\tffluid_arrayfunc_module ffluid_data_copy ffluid_data_fma\n");
}

void ffluid_data_copy(data_ptr in, data_ptr out) {
  memcpy(out->Q, in->Q, in->N*sizeof(long_complex_t));
  memcpy(out->V, in->V, in->N*sizeof(long_complex_t));
  out->q0 = in->q0;
  out->u0 = in->u0;
  out->l = in->l;
  out->time = in->time;
}

void ffluid_data_init_copy(data_ptr in, data_ptr out) {
  out->N = in->N;
  ffluid_init_data(out);
  ffluid_data_copy(in, out);
  ffluid_setup_grid(out);
}


void ffluid_data_fma(double op1, data_ptr op2, data_ptr op3, data_ptr out) {
}
